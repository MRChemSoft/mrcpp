#include "Printer.h"
#include "Timer.h"

#include "Bank.h"

namespace mrcpp {

using namespace Eigen;
using namespace std;

#ifdef MRCPP_HAS_MPI
    int metadata_block[3]; // can add more metadata in future
    int const size_metadata = 3;
#endif

Bank::~Bank() {
    // delete all data and accounts
}

struct Blockdata_struct {
    std::vector<double *> data; // to store the incoming data. One column for each orbital on the same node.
    int N_rows = 0;             // the number of coefficients in one column of the block.
    std::map<int, int> id2data; // internal index of the data in the block
    std::vector<int> id;        // the id of each column. Either nodeid, or orbid
};
struct OrbBlock_struct {
    std::vector<double *> data; // pointer to the data
    std::map<int, int> id2data; // internal index of the data in the block
    std::vector<int> id;        // the nodeid of the data
    // note that N_rows can be different inside the same orbblock: root node have scaling and wavelets, other nodes have only wavelets
};
struct mem_struct {
    std::vector<double *> chunk_p; // vector with allocated chunks
    int p = -1;                    // position of next available memory (not allocated if < 0)
    // on Betzy 1024*1024*4 ok, 1024*1024*2 NOT ok: leads to memory fragmentation (on "Betzy" 2023)
    int chunk_size = 1024 * 1024 * 4; // chunksize (in number of doubles). data_p[i]+chunk_size is end of chunk i
    int account = -1;
    double *get_mem(int size) {
        if (p < 0 or size > chunk_size or p + size > chunk_size) { // allocate new chunk of memory
            if (size > 1024 * 1024) {
                // make a special chunk just for this
                double *m_p = new double[size];
                chunk_p.push_back(m_p);
                p = -1;
                return m_p;
            } else {
                double *m_p = new double[chunk_size];
                chunk_p.push_back(m_p);
                p = 0;
            }
        }
        double *m_p = chunk_p[chunk_p.size() - 1] + p;
        p += size;
        return m_p;
    }
};
std::map<int, std::map<int, Blockdata_struct> *> get_nodeid2block; // to get block from its nodeid (all coeff for one node)
std::map<int, std::map<int, OrbBlock_struct> *> get_orbid2block;   // to get block from its orbid

std::map<int, mem_struct *> mem;

#ifdef MRCPP_HAS_MPI
int const MIN_SCALE = -999; // Smaller than smallest scale
#endif
int naccounts = 0;

void Bank::open() {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    char safe_data1;
    int deposit_size = sizeof(deposit);
    int n_chunks, ix;
    int messages[message_size];
    int datasize = -1;
    std::map<int, int> get_numberofclients;
    std::map<NodeIndex<3>, int> nIdx2id;
    bool printinfo = false;
    int max_account_id = -1;
    int next_task = 0;
    int tot_ntasks = 0;
    std::map<int, std::vector<int>> readytasks;
    // The bank never goes out of this loop until it receives a close message!
    while (true) {
        MPI_Recv(messages, message_size, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_bank, &status);
        if (printinfo)
            std::cout << world_rank << " got message " << messages[0] << " from " << status.MPI_SOURCE << " account " << messages[1] << " last account " << max_account_id << " message2 "
                      << messages[2] << std::endl;
        int message = messages[0];

        // can be called directly:
        if (message == CLOSE_BANK) {
            if (is_bank and printinfo) std::cout << "Bank is closing" << std::endl;
            this->clear_bank();
            break; // close bank, i.e stop listening for incoming messages
        } else if (message == GET_MAXTOTDATA) {
            int maxsize_int = maxsize / 1024; // convert into MB
            MPI_Send(&maxsize_int, 1, MPI_INT, status.MPI_SOURCE, 1171, comm_bank);
            continue;
        } else if (message == GET_TOTDATA) {
            int maxsize_int = totcurrentsize / 1024; // convert into MB
            MPI_Send(&maxsize_int, 1, MPI_INT, status.MPI_SOURCE, 1172, comm_bank);
            continue;
        } else if (message == NEW_ACCOUNT) {
            // we just have to pick out a number that is not already assigned
            int account = (max_account_id + 1) % 1000000000;
            while (get_deposits.count(account)) account = (account + 1) % 1000000000; // improbable this is used
            max_account_id = account;
            naccounts++;
            // create default content
            get_deposits[account] = new std::vector<deposit>;
            get_deposits[account]->resize(1);
            get_id2ix[account] = new std::map<int, int>;
            get_id2qu[account] = new std::map<int, int>;
            get_queue[account] = new std::vector<queue_struct>;
            get_queue[account]->resize(1);
            get_orbid2block[account] = new std::map<int, OrbBlock_struct>;
            get_nodeid2block[account] = new std::map<int, Blockdata_struct>;
            get_numberofclients[account] = messages[1];
            get_readytasks[account] = new std::map<int, std::vector<int>>;
            currentsize[account] = 0;
            mem[account] = new mem_struct;
            mem[account]->account = account;
            MPI_Send(&account, 1, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
            continue;
        }

        // the following is only accessible through an account

        int account = messages[1];
        auto it_dep = get_deposits.find(account);
        if (it_dep == get_deposits.end() || it_dep->second == nullptr) {
            cout << "ERROR, my dep account does not exist!! " << account << " " << message << endl;
            MSG_ABORT("Account error");
        }
        std::vector<deposit> &deposits = *get_deposits[account];
        std::map<int, int> &id2ix = *get_id2ix[account]; // gives zero if id is not defined
        std::map<int, int> &id2qu = *get_id2qu[account];
        std::vector<queue_struct> &queue = *get_queue[account];
        std::map<int, OrbBlock_struct> &orbid2block = *get_orbid2block[account];
        std::map<int, Blockdata_struct> &nodeid2block = *get_nodeid2block[account];
        auto it_tasks = get_readytasks.find(account);
        if (it_tasks == get_readytasks.end() || it_tasks->second == nullptr) {
            cout << "ERROR, my account does not exist!! " << account << " " << message << endl;
            MSG_ABORT("Account error");
        }
        std::map<int, std::vector<int>> &readytasks = *get_readytasks[account];

        if (message == CLOSE_ACCOUNT) {
            get_numberofclients[account]--;
            if (get_numberofclients[account] == 0) {
                // all clients have closed the account. We remove the account.
                remove_account(account);
            }
        }

        else if (message == CLEAR_BANK) {
            this->clear_bank();
            for (auto const &block : nodeid2block) {
                if (block.second.data.size() > 0) {
                    currentsize[account] -= block.second.N_rows * block.second.data.size() / 128; // converted into kB
                    totcurrentsize -= block.second.N_rows * block.second.data.size() / 128;       // converted into kB
                }
            }
            nodeid2block.clear();
            orbid2block.clear();
            // send message that it is ready (value of message is not used)
            MPI_Ssend(&message, 1, MPI_INT, status.MPI_SOURCE, 77, comm_bank);
        }

        else if (message == GET_NODEDATA or message == GET_NODEBLOCK) {
            // NB: has no queue system yet
            int nodeid = messages[2]; // which block to fetch from
            if (nodeid2block.count(nodeid)) {
                Blockdata_struct &block = nodeid2block[nodeid];
                int dataindex = 0; // internal index of the data in the block
                int size = 0;
                if (message == GET_NODEDATA) {
                    int orbid = messages[3];          // which part of the block to fetch
                    dataindex = block.id2data[orbid]; // column of the data in the block
                    size = block.N_rows;              // number of doubles to fetch
                    if (size != messages[4]) std::cout << "ERROR nodedata has wrong size" << std::endl;
                    double *data_p = block.data[dataindex];
                    if (size > 0) MPI_Send(data_p, size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
                } else {
                    // send entire block. First make one contiguous superblock
                    // Prepare the data as one contiguous block
                    if (block.data.size() == 0) std::cout << "Zero size blockdata! " << nodeid << " " << block.N_rows << std::endl;
                    MatrixXd DataBlock(block.N_rows, block.data.size());
                    size = block.N_rows * block.data.size();
                    if (printinfo) std::cout << " rewrite into superblock " << block.data.size() << " " << block.N_rows << " nodeid " << nodeid << std::endl;
                    for (int j = 0; j < block.data.size(); j++) {
                        for (int i = 0; i < block.N_rows; i++) { DataBlock(i, j) = block.data[j][i]; }
                    }
                    dataindex = 0; // start from first column
                    // send info about the size of the superblock
                    metadata_block[0] = nodeid;            // nodeid
                    metadata_block[1] = block.data.size(); // number of columns
                    metadata_block[2] = size;              // total size = rows*columns
                    MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
                    // send info about the id of each column
                    MPI_Send(block.id.data(), metadata_block[1], MPI_INT, status.MPI_SOURCE, 2, comm_bank);
                    if (size > 0) MPI_Send(DataBlock.data(), size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
                }
            } else {
                if (printinfo) std::cout << " block " << nodeid << " does not exist " << std::endl;
                // Block with this id does not exist.
                if (message == GET_NODEDATA) {
                    int size = messages[4]; // number of doubles to send
                    if (size == 0) {
                        std::cout << "WARNING: GET_NODEDATA asks for zero size data" << std::endl;
                        metadata_block[2] = size;
                        MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 3, comm_bank);
                    } else {
                        std::vector<double> zero(size, 0.0); // send zeroes
                        MPI_Ssend(zero.data(), size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
                    }
                } else {
                    metadata_block[0] = nodeid;
                    metadata_block[1] = 0; // number of columns
                    metadata_block[2] = 0; // total size = rows*columns
                    MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
                }
            }
        } else if (message == GET_ORBBLOCK) {
            // NB: BLOCKDATA has no queue system yet
            int orbid = messages[2]; // which block to fetch from

            if (orbid2block.count(orbid)) {
                OrbBlock_struct &block = orbid2block[orbid];
                if (block.data.size() == 0) std::cout << "Zero size blockdata! C " << orbid << " " << std::endl;
                // send entire block. First make one contiguous superblock
                // Prepare the data as one contiguous block
                int size = 0;
                for (int j = 0; j < block.data.size(); j++) {
                    int nodeid = block.id[j];
                    int Nrows = nodeid2block[nodeid].N_rows; // note that root nodes have scaling and wavelets, while other nodes have only wavelets -> N_rows is not a constant.
                    size += Nrows;
                }
                std::vector<double> coeff(size);
                int ij = 0;
                for (int j = 0; j < block.data.size(); j++) {
                    int nodeid = block.id[j];
                    int Nrows = nodeid2block[nodeid].N_rows;
                    for (int i = 0; i < Nrows; i++) { coeff[ij++] = block.data[j][i]; }
                }
                // send info about the size of the superblock
                metadata_block[0] = orbid;
                metadata_block[1] = block.data.size(); // number of columns
                metadata_block[2] = size;              // total size = rows*columns
                MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
                MPI_Send(block.id.data(), metadata_block[1], MPI_INT, status.MPI_SOURCE, 2, comm_bank);
                MPI_Send(coeff.data(), size, MPI_DOUBLE, status.MPI_SOURCE, 3, comm_bank);
            } else {
                // it is possible and allowed that the block has not been written
                if (printinfo) std::cout << " block does not exist " << orbid << " " << orbid2block.count(orbid) << std::endl;
                // Block with this id does not exist.
                metadata_block[0] = orbid;
                metadata_block[1] = 0; // number of columns
                metadata_block[2] = 0; // total size = rows*columns
                MPI_Send(metadata_block, size_metadata, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
            }
        }

        else if (message == GET_FUNCTION or message == GET_FUNCTION_AND_WAIT or message == GET_FUNCTION_AND_DELETE or message == GET_FUNCTION or message == GET_DATA) {
            // withdrawal
            int id = messages[2];
            if (message == GET_DATA and messages[3] > MIN_SCALE) {
                NodeIndex<3> nIdx;
                nIdx.setScale(messages[4]);
                nIdx.setTranslation({messages[2], messages[5], messages[6]});
                if (nIdx2id.count(nIdx) == 0) {
                    // data is not yet saved, but one can hope it will be created at some stage
                    id = nIdx2id.size();
                    nIdx2id[nIdx] = id;
                } else {
                    id = nIdx2id[nIdx];
                }
            }
            int ix = id2ix[id];
            if (id2ix.count(id) == 0 or ix == 0) {
                if (printinfo) std::cout << world_rank << " not found " << id << " " << message << std::endl;
                if (message == GET_FUNCTION or message == GET_FUNCTION_AND_DELETE) {
                    // do not wait for the orbital to arrive
                    int found = 0;
                    if (printinfo) std::cout << world_rank << " sending found 0 to " << status.MPI_SOURCE << std::endl;
                    MPI_Send(&found, 1, MPI_INT, status.MPI_SOURCE, 117, comm_bank);
                } else {
                    // the id does not exist. Put in queue and Wait until it is defined
                    if (printinfo) std::cout << world_rank << " queuing " << id << " " << id2ix.count(id) << ix << std::endl;
                    if (id2qu[id] == 0) {
                        queue.push_back({id, {status.MPI_SOURCE}});
                        id2qu[id] = queue.size() - 1;
                    } else {
                        // somebody is already waiting for this id. queue in queue
                        queue[id2qu[id]].clients.push_back(status.MPI_SOURCE);
                    }
                }
            } else {
                int ix = id2ix[id];
                if (deposits[ix].id != id) std::cout << ix << " Bank accounting error " << std::endl;
                if (message == GET_FUNCTION or message == GET_FUNCTION_AND_WAIT or message == GET_FUNCTION_AND_DELETE) {
                    if (message == GET_FUNCTION or message == GET_FUNCTION_AND_DELETE) {
                        int found = 1;
                        MPI_Send(&found, 1, MPI_INT, status.MPI_SOURCE, 117, comm_bank);
                    }
                    send_function(*deposits[ix].orb, status.MPI_SOURCE, 1, comm_bank);
                    if (message == GET_FUNCTION_AND_DELETE) {
                        currentsize[account] -= deposits[ix].orb->getSizeNodes();
                        totcurrentsize -= deposits[ix].orb->getSizeNodes();
                        deposits[ix].orb->free();
                        id2ix[id] = 0;
                    }
                }
                if (message == GET_DATA) { MPI_Send(deposits[ix].data, deposits[ix].datasize, MPI_DOUBLE, status.MPI_SOURCE, 1, comm_bank); }
            }
        } else if (message == SAVE_NODEDATA) {
            int nodeid = messages[2]; // which block to write
            int orbid = messages[3];  // which part of the block
            int size = messages[4];   // number of doubles

            // test if the block exists already
            if (printinfo) std::cout << world_rank << " save data nodeid " << nodeid << " size " << size << std::endl;
            // append the incoming data
            Blockdata_struct &block = nodeid2block[nodeid];
            block.id2data[orbid] = nodeid2block[nodeid].data.size(); // internal index of the data in the block
            double *data_p = mem[account]->get_mem(size);            // new double[size];
            currentsize[account] += size / 128;                      // converted into kB
            totcurrentsize += size / 128;                            // converted into kB
            this->maxsize = std::max(totcurrentsize, this->maxsize);
            block.data.push_back(data_p);
            block.id.push_back(orbid);
            if (block.N_rows > 0 and block.N_rows != size) cout << " ERROR block size incompatible " << block.N_rows << " " << size << endl;
            block.N_rows = size;

            OrbBlock_struct &orbblock = orbid2block[orbid];
            orbblock.id2data[nodeid] = orbblock.data.size(); // internal index of the data in the block
            orbblock.data.push_back(data_p);
            orbblock.id.push_back(nodeid);
            // orbblock.N_rows.push_back(size);

            MPI_Recv(data_p, size, MPI_DOUBLE, status.MPI_SOURCE, 1, comm_bank, &status);
            if (printinfo) std::cout << " written block " << nodeid << " id " << orbid << " subblocks " << nodeid2block[nodeid].data.size() << std::endl;
        } else if (message == SAVE_FUNCTION or message == SAVE_DATA) {
            // make a new deposit
            int id = messages[2];
            if (message == SAVE_DATA and messages[4] > MIN_SCALE) {
                // has to find or create unique id from NodeIndex. Use the same internal mapping for all trees
                NodeIndex<3> nIdx;
                nIdx.setScale(messages[4]);
                nIdx.setTranslation({messages[2], messages[5], messages[6]});
                if (nIdx2id.count(nIdx) == 0) {
                    id = nIdx2id.size();
                    nIdx2id[nIdx] = id;
                } else {
                    id = nIdx2id[nIdx];
                }
            }
            int exist_flag = 0;
            if (id2ix[id]) {
                std::cout << "WARNING: id " << id << " exists already"
                          << " " << status.MPI_SOURCE << " " << message << " " << messages[1] << std::endl;
                ix = id2ix[id]; // the deposit exist from before. Will be overwritten
                exist_flag = 1;
                if (message == SAVE_DATA and !deposits[ix].hasdata) {
                    datasize = messages[3];
                    exist_flag = 0;
                    //                    deposits[ix].data = new double[datasize];
                    deposits[ix].data = mem[account]->get_mem(datasize);
                    currentsize[account] += datasize / 128; // converted into kB
                    totcurrentsize += datasize / 128;       // converted into kB
                    this->maxsize = std::max(totcurrentsize, this->maxsize);
                    deposits[ix].hasdata = true;
                }
            } else {
                ix = deposits.size(); // NB: ix is now index of last element + 1
                deposits.resize(ix + 1);
                if (message == SAVE_FUNCTION) deposits[ix].orb = new CompFunction<3>(0);
                if (message == SAVE_DATA) {
                    datasize = messages[3];
                    deposits[ix].data = mem[account]->get_mem(datasize); // new double[datasize];
                    currentsize[account] += datasize / 128;              // converted into kB
                    totcurrentsize += datasize / 128;                    // converted into kB
                    this->maxsize = std::max(totcurrentsize, this->maxsize);
                    deposits[ix].hasdata = true;
                }
            }
            deposits[ix].id = id;
            id2ix[deposits[ix].id] = ix;
            deposits[ix].source = status.MPI_SOURCE;
            if (message == SAVE_FUNCTION) {
                recv_function(*deposits[ix].orb, deposits[ix].source, 1, comm_bank);
                if (exist_flag == 0) {
                    currentsize[account] += deposits[ix].orb->getSizeNodes();
                    totcurrentsize += deposits[ix].orb->getSizeNodes();
                    this->maxsize = std::max(totcurrentsize, this->maxsize);
                }
            }
            if (message == SAVE_DATA) {
                datasize = messages[3];
                deposits[ix].datasize = datasize;
                MPI_Recv(deposits[ix].data, datasize, MPI_DOUBLE, deposits[ix].source, 1, comm_bank, &status);
            }
            if (id2qu[deposits[ix].id] != 0) {
                // someone is waiting for those data. Send to them
                int iq = id2qu[deposits[ix].id];
                if (deposits[ix].id != queue[iq].id) std::cout << ix << " Bank queue accounting error " << std::endl;
                for (int iqq : queue[iq].clients) {
                    if (message == SAVE_FUNCTION) { send_function(*deposits[ix].orb, iqq, 1, comm_bank); }
                    if (message == SAVE_DATA) { MPI_Send(deposits[ix].data, messages[3], MPI_DOUBLE, iqq, 1, comm_bank); }
                }
                queue[iq].clients.clear(); // cannot erase entire queue[iq], because that would require to shift all the
                                           // id2qu value larger than iq
                queue[iq].id = -1;
                id2qu.erase(deposits[ix].id);
            }

            // Task manager members:
        } else if (message == INIT_TASKS) {
            tot_ntasks = messages[2];
            next_task = 0;
        } else if (message == GET_NEXTTASK) {
            int task = next_task;
            if (next_task >= tot_ntasks) task = -1; // flag to show all tasks are assigned
            MPI_Send(&task, 1, MPI_INT, status.MPI_SOURCE, 1, comm_bank);
            next_task++;
        } else if (message == PUT_READYTASK) {
            readytasks[messages[2]].push_back(messages[3]);
        }
        if (message == DEL_READYTASK) {
            for (int i = 0; i < readytasks[messages[2]].size(); i++) { // we expect small sizes
                if (readytasks[messages[2]][i] == messages[3]) {
                    readytasks[messages[2]].erase(readytasks[messages[2]].begin() + i);
                    break;
                }
            }
        } else if (message == GET_READYTASK) {
            int nready = 0;
            if (readytasks.count(messages[2]) > 0) nready = readytasks[messages[2]].size();
            MPI_Send(&nready, 1, MPI_INT, status.MPI_SOURCE, 844, mpi::comm_bank);
            if (nready > 0) MPI_Send(readytasks[messages[2]].data(), nready, MPI_INT, status.MPI_SOURCE, 845, mpi::comm_bank);
        } else if (message == GET_READYTASK_DEL) {
            int nready = 0;
            if (readytasks.count(messages[2]) > 0) { nready = readytasks[messages[2]].size(); }
            MPI_Send(&nready, 1, MPI_INT, status.MPI_SOURCE, 844, mpi::comm_bank);
            if (nready > 0) MPI_Send(readytasks[messages[2]].data(), nready, MPI_INT, status.MPI_SOURCE, 845, mpi::comm_bank);
            if (nready > 0) readytasks[messages[2]].resize(0);
        }
    }
#endif
}

// Ask to close the Bank
void Bank::close() {
#ifdef MRCPP_HAS_MPI
    int messages[message_size];
    messages[0] = CLOSE_BANK;
    for (int i = 0; i < bank_size; i++) MPI_Send(messages, 1, MPI_INT, bankmaster[i], 0, comm_bank);
    if (tot_bank_size > bank_size) MPI_Send(messages, 1, MPI_INT, task_bank, 0, comm_bank);
#endif
}

void Bank::clear_bank() {
#ifdef MRCPP_HAS_MPI
    for (auto account : accounts) { remove_account(account); }
#endif
}

int Bank::clearAccount(int account, int iclient, MPI_Comm comm) {
#ifdef MRCPP_HAS_MPI
    closeAccount(account);
    return openAccount(iclient, comm);
#else
    return 1;
#endif
}
void Bank::remove_account(int account) {
#ifdef MRCPP_HAS_MPI
    naccounts--;
    auto it = get_deposits.find(account);
    if (it == get_deposits.end() || it->second == nullptr) {
        cout << "ERROR, my account depositsdoes not exist!! " << account << " " << endl;
        MSG_ABORT("depositsAccount error");
    }
    std::vector<deposit> &deposits = *get_deposits[account];
    for (int ix = 1; ix < deposits.size(); ix++) {
        if (deposits[ix].orb != nullptr) deposits[ix].orb->free();
        if (deposits[ix].hasdata) {
            currentsize[account] -= deposits[ix].datasize / 128;
            totcurrentsize -= deposits[ix].datasize / 128;
        }
        if (deposits[ix].hasdata) (*get_id2ix[account])[deposits[ix].id] = 0; // indicate that it does not exist
        deposits[ix].hasdata = false;
    }
    deposits.clear();
    get_deposits.erase(account);
    delete get_queue[account];
    get_queue.erase(account);
    delete get_id2ix[account];
    get_id2ix.erase(account);
    delete get_id2qu[account];
    get_id2qu.erase(account);
    delete get_readytasks[account];
    get_readytasks.erase(account);

    std::map<int, Blockdata_struct> &nodeid2block = *get_nodeid2block[account];
    std::map<int, OrbBlock_struct> &orbid2block = *get_orbid2block[account];

    for (auto const &block : nodeid2block) {
        currentsize[account] -= block.second.N_rows * block.second.data.size() / 128; // converted into kB
        totcurrentsize -= block.second.N_rows * block.second.data.size() / 128;       // converted into kB
    }
    nodeid2block.clear();
    orbid2block.clear();

    get_nodeid2block.erase(account);
    get_orbid2block.erase(account);

    for (double *c_p : mem[account]->chunk_p) delete[] c_p;
    mem.erase(account);
    currentsize.erase(account);
#endif
}

int Bank::openAccount(int iclient, MPI_Comm comm) {
    // NB: this is a collective call, since we need all the accounts to be synchronized
    int account_id[1] = {-1};
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = NEW_ACCOUNT;
    messages[1] = 0;
    int size;
    MPI_Comm_size(comm, &size);
    messages[1] = size;
    if (iclient == 0) {
        for (int i = 0; i < bank_size; i++) {
            int account_id_i[1];
            MPI_Send(messages, message_size, MPI_INT, bankmaster[i], 0, comm_bank);
            MPI_Recv(account_id_i, 1, MPI_INT, bankmaster[i], 1, comm_bank, &status);
            if (i > 0 and account_id_i[0] != account_id[0]) MSG_ABORT("Account id mismatch!");
            account_id[0] = account_id_i[0];
        }
        MPI_Bcast(account_id, 1, MPI_INT, 0, comm);
    } else {
        MPI_Bcast(account_id, 1, MPI_INT, 0, comm);
    }
#endif
    return account_id[0];
}

int Bank::openTaskManager(int ntasks, int iclient, MPI_Comm comm) {
    // NB: this is a collective call, since we need all the accounts to be synchronized
    int account_id = -1;
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = NEW_ACCOUNT;
    int size;
    MPI_Comm_size(comm, &size);
    messages[1] = size;
    if (iclient == 0) {
        MPI_Send(messages, 2, MPI_INT, task_bank, 0, comm_bank);
        MPI_Recv(&account_id, 1, MPI_INT, task_bank, 1, comm_bank, &status);
        if (tot_bank_size == bank_size) {
            // make a dummy account so that all account_id are synchronized
            int account_id_i;
            for (int i = 0; i < bank_size; i++) {
                if (bankmaster[i] != task_bank) {
                    MPI_Send(messages, 2, MPI_INT, bankmaster[i], 0, comm_bank);
                    MPI_Recv(&account_id_i, 1, MPI_INT, bankmaster[i], 1, comm_bank, &status);
                    if (i > 0 and account_id_i != account_id) MSG_ABORT("Account id mismatch!");
                }
            }
        }
        messages[0] = INIT_TASKS;
        messages[1] = account_id;
        messages[2] = ntasks;
        MPI_Send(messages, 3, MPI_INT, task_bank, 2, comm_bank);
        MPI_Bcast(&account_id, 1, MPI_INT, 0, comm);
    } else {
        MPI_Bcast(&account_id, 1, MPI_INT, 0, comm);
    }
#endif
    return account_id;
}

void Bank::closeAccount(int account_id) {
// The account will in reality not be removed before everybody has sent a close message
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = CLOSE_ACCOUNT;
    messages[1] = account_id;
    for (int i = 0; i < bank_size; i++) MPI_Send(messages, 2, MPI_INT, bankmaster[i], 0, comm_bank);
#endif
}

void Bank::closeTaskManager(int account_id) {
// The account will in reality not be removed before everybody has sent a close message
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = CLOSE_ACCOUNT;
    messages[1] = account_id;
    MPI_Send(messages, 2, MPI_INT, task_bank, 0, comm_bank);
#endif
}

int Bank::get_maxtotalsize() {
    int maxtot = 0;
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int datasize;
    int messages[message_size];
    messages[0] = GET_MAXTOTDATA;
    for (int i = 0; i < bank_size; i++) {
        MPI_Send(messages, 1, MPI_INT, bankmaster[i], 0, comm_bank);
        MPI_Recv(&datasize, 1, MPI_INT, bankmaster[i], 1171, comm_bank, &status);
        maxtot = std::max(maxtot, datasize);
    }
#endif
    return maxtot;
}

std::vector<int> Bank::get_totalsize() {
    std::vector<int> tot;
#ifdef HAVE_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_TOTDATA;
    int datasize;
    for (int i = 0; i < bank_size; i++) {
        MPI_Send(messages, 1, MPI_INT, bankmaster[i], 0, comm_bank);
        MPI_Recv(&datasize, 1, MPI_INT, bankmaster[i], 1172, comm_bank, &status);
        tot.push_back(datasize);
    }
#endif
    return tot;
}

// Accounts: (clients)

// save orbital in Bank with identity id

// get orbital with identity id.
// If wait=0, return immediately with value zero if not available (default)
// else, wait until available
int BankAccount::get_func(int id, CompFunction<3> &func, int wait) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[1] = account_id;
    messages[2] = id;
    if (wait == 0) {
        messages[0] = GET_FUNCTION;
        MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
        int found;
        MPI_Recv(&found, 1, MPI_INT, bankmaster[id % bank_size], 117, comm_bank, &status);
        if (found != 0) {
            recv_function(func, bankmaster[id % bank_size], 1, comm_bank);
            return 1;
        } else {
            return 0;
        }
    } else {
        messages[0] = GET_FUNCTION_AND_WAIT;
        MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
        recv_function(func, bankmaster[id % bank_size], 1, comm_bank);
    }
#endif
    return 1;
}

// get orbital with identity id, and delete from bank.
// return immediately with value zero if not available
int BankAccount::get_func_del(int id, CompFunction<3> &orb) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_FUNCTION_AND_DELETE;
    messages[1] = account_id;
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    int found;
    MPI_Recv(&found, 1, MPI_INT, bankmaster[id % bank_size], 117, comm_bank, &status);
    if (found != 0) {
        recv_function(orb, bankmaster[id % bank_size], 1, comm_bank);
        return 1;
    } else {
        return 0;
    }
#endif
    return 1;
}

// save function in Bank with identity id
int BankAccount::put_func(int id, CompFunction<3> &func) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to id
    int messages[message_size];
    messages[0] = SAVE_FUNCTION;
    messages[1] = account_id;
    messages[2] = id;
    MPI_Send(messages, 3, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    send_function(func, bankmaster[id % bank_size], 1, comm_bank);
#endif
    return 1;
}

// save data in Bank with identity id . datasize MUST have been set already. NB:not tested
int BankAccount::put_data(int id, int size, double *data) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to id
    int messages[message_size];

    messages[0] = SAVE_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = size;
    messages[4] = MIN_SCALE; // to indicate that it is defined by id
    MPI_Send(messages, 5, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank);
#endif
    return 1;
}

// save data in Bank with identity id . datasize MUST have been set already. NB:not tested
int BankAccount::put_data(int id, int size, ComplexDouble *data) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to id
    int messages[message_size];

    messages[0] = SAVE_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = size * 2;  // save as twice as many doubles
    messages[4] = MIN_SCALE; // to indicate that it is defined by id
    MPI_Send(messages, 5, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank);
#endif
    return 1;
}

// save data in Bank with identity nIdx. datasize MUST have been set already. NB:not tested
int BankAccount::put_data(NodeIndex<3> nIdx, int size, double *data) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to id
    int messages[message_size];
    messages[0] = SAVE_DATA;
    messages[1] = account_id;
    messages[2] = nIdx.getTranslation(0);
    messages[3] = size;
    messages[4] = nIdx.getScale();
    messages[5] = nIdx.getTranslation(1);
    messages[6] = nIdx.getTranslation(2);
    int id = std::abs(nIdx.getTranslation(0) + nIdx.getTranslation(1) + nIdx.getTranslation(2));
    MPI_Send(messages, 7, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank);
#endif
    return 1;
}

// save data in Bank with identity nIdx. datasize MUST have been set already. NB:not tested
int BankAccount::put_data(NodeIndex<3> nIdx, int size, ComplexDouble *data) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to id
    int messages[message_size];
    messages[0] = SAVE_DATA;
    messages[1] = account_id;
    messages[2] = nIdx.getTranslation(0);
    messages[3] = size * 2; // save as twice as many doubles
    messages[4] = nIdx.getScale();
    messages[5] = nIdx.getTranslation(1);
    messages[6] = nIdx.getTranslation(2);
    int id = std::abs(nIdx.getTranslation(0) + nIdx.getTranslation(1) + nIdx.getTranslation(2));
    MPI_Send(messages, 7, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_data(int id, int size, double *data) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = MIN_SCALE;
    MPI_Send(messages, 4, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank, &status);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_data(int id, int size, ComplexDouble *data) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = MIN_SCALE;
    MPI_Send(messages, 4, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    // fetch as twice as many doubles
    MPI_Recv(data, size * 2, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank, &status);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_data(NodeIndex<3> nIdx, int size, double *data) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    int id = std::abs(nIdx.getTranslation(0) + nIdx.getTranslation(1) + nIdx.getTranslation(2));
    messages[0] = GET_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = nIdx.getScale();
    messages[4] = nIdx.getTranslation(0);
    messages[5] = nIdx.getTranslation(1);
    messages[6] = nIdx.getTranslation(2);
    MPI_Send(messages, 7, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank, &status);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_data(NodeIndex<3> nIdx, int size, ComplexDouble *data) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int messages[message_size];
    int id = std::abs(nIdx.getTranslation(0) + nIdx.getTranslation(1) + nIdx.getTranslation(2));
    messages[0] = GET_DATA;
    messages[1] = account_id;
    messages[2] = id;
    messages[3] = nIdx.getScale();
    messages[4] = nIdx.getTranslation(0);
    messages[5] = nIdx.getTranslation(1);
    messages[6] = nIdx.getTranslation(2);
    MPI_Send(messages, 7, MPI_INT, bankmaster[id % bank_size], 0, comm_bank);
    // fetch as twice as many doubles
    MPI_Recv(data, size * 2, MPI_DOUBLE, bankmaster[id % bank_size], 1, comm_bank, &status);
#endif
    return 1;
}

// save data in Bank with identity id as part of block with identity nodeid.
int BankAccount::put_nodedata(int id, int nodeid, int size, double *data) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to nodeid
    int messages[message_size];
    messages[0] = SAVE_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block
    messages[4] = size;   // size of this data
    MPI_Send(messages, 5, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Send(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 1, comm_bank);
#endif
    return 1;
}

// save data in Bank with identity id as part of block with identity nodeid.
// NB: Complex is stored as two doubles
int BankAccount::put_nodedata(int id, int nodeid, int size, ComplexDouble *data) {
#ifdef MRCPP_HAS_MPI
    // for now we distribute according to nodeid
    int messages[message_size];
    messages[0] = SAVE_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid;   // which block
    messages[3] = id;       // id within block
    messages[4] = 2 * size; // size of this data
    MPI_Send(messages, 5, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Send(data, 2 * size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 1, comm_bank);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_nodedata(int id, int nodeid, int size, double *data, std::vector<int> &idVec) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    // get the column with identity id
    int messages[message_size];
    messages[0] = GET_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block.
    messages[4] = size;   // expected size of data
    MPI_Send(messages, 5, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
#endif
    return 1;
}

// get data with identity id
int BankAccount::get_nodedata(int id, int nodeid, int size, ComplexDouble *data, std::vector<int> &idVec) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    // get the column with identity id
    int messages[message_size];
    messages[0] = GET_NODEDATA;
    messages[1] = account_id;
    messages[2] = nodeid; // which block
    messages[3] = id;     // id within block.
    messages[4] = size;   // expected size of data
    MPI_Send(messages, 5, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
#endif
    return 1;
}

// get all data for nodeid (same nodeid, different orbitals)
int BankAccount::get_nodeblock(int nodeid, double *data, std::vector<int> &idVec) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    // get the entire superblock and also the id of each column
    int messages[message_size];
    messages[0] = GET_NODEBLOCK;
    messages[1] = account_id;
    messages[2] = nodeid;

    MPI_Send(messages, 3, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], 1, comm_bank, &status);
    idVec.resize(metadata_block[1]);
    int size = metadata_block[2];
    if (size > 0) MPI_Recv(idVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], 2, comm_bank, &status);
    if (size > 0) MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
#endif
    return 1;
}

// get all data for nodeid (same nodeid, different orbitals)
int BankAccount::get_nodeblock(int nodeid, ComplexDouble *data, std::vector<int> &idVec) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    // get the entire superblock and also the id of each column
    int messages[message_size];
    messages[0] = GET_NODEBLOCK;
    messages[1] = account_id;
    messages[2] = nodeid;

    MPI_Send(messages, 3, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], 1, comm_bank, &status);
    idVec.resize(metadata_block[1]);
    int size = metadata_block[2];
    if (size > 0) MPI_Recv(idVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], 2, comm_bank, &status);
    if (size > 0) MPI_Recv(data, size, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
#endif
    return 1;
}

// get all data with identity orbid (same orbital, different nodes)
int BankAccount::get_orbblock(int orbid, double *&data, std::vector<int> &nodeidVec, int bankstart) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int nodeid = wrk_rank + bankstart;
    // get the entire superblock and also the nodeid of each column
    int messages[message_size];
    messages[0] = GET_ORBBLOCK;
    messages[1] = account_id;
    messages[2] = orbid;
    MPI_Send(messages, 3, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], 1, comm_bank, &status);
    nodeidVec.resize(metadata_block[1]);
    int totsize = metadata_block[2];
    if (totsize > 0) MPI_Recv(nodeidVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], 2, comm_bank, &status);
    data = new double[totsize];
    if (totsize > 0) MPI_Recv(data, totsize, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
#endif
    return 1;
}

// get all data with identity orbid (same orbital, different nodes)
int BankAccount::get_orbblock(int orbid, ComplexDouble *&data, std::vector<int> &nodeidVec, int bankstart) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    int nodeid = wrk_rank + bankstart;
    // get the entire superblock and also the nodeid of each column
    int messages[message_size];
    messages[0] = GET_ORBBLOCK;
    messages[1] = account_id;
    messages[2] = orbid;
    MPI_Send(messages, 3, MPI_INT, bankmaster[nodeid % bank_size], 0, comm_bank);
    MPI_Recv(metadata_block, size_metadata, MPI_INT, bankmaster[nodeid % bank_size], 1, comm_bank, &status);
    nodeidVec.resize(metadata_block[1]);
    int totsize = metadata_block[2];
    if (totsize > 0) MPI_Recv(nodeidVec.data(), metadata_block[1], MPI_INT, bankmaster[nodeid % bank_size], 2, comm_bank, &status);
    data = new ComplexDouble[totsize / 2];
    if (totsize > 0) MPI_Recv(data, totsize, MPI_DOUBLE, bankmaster[nodeid % bank_size], 3, comm_bank, &status);
#endif
    return 1;
}

// creator. NB: collective
BankAccount::BankAccount(int iclient, MPI_Comm comm) {
    this->account_id = dataBank.openAccount(iclient, comm);
#ifdef MRCPP_HAS_MPI
    MPI_Barrier(comm);
#endif
}

// destructor
BankAccount::~BankAccount() {
    // The account will in reality not be removed before everybody has sent a delete message
    dataBank.closeAccount(this->account_id);
}

// closes account and reopen a new empty account. NB: account_id will change
void BankAccount::clear(int iclient, MPI_Comm comm) {
    this->account_id = dataBank.clearAccount(this->account_id, iclient, comm);
}

// creator. NB: collective
TaskManager::TaskManager(int ntasks, int iclient, MPI_Comm comm) {
    this->n_tasks = ntasks;
    if (bank_size == 0) return;
    this->account_id = dataBank.openTaskManager(ntasks, iclient, comm);
#ifdef MRCPP_HAS_MPI
    MPI_Barrier(comm);
#endif
}

// destructor
TaskManager::~TaskManager() {
    // The account will in reality not be removed before everybody has sent a delete message
    if (this->account_id < 0) return;
    dataBank.closeTaskManager(this->account_id);
}

int TaskManager::next_task() {
    int nexttask = 0;
#ifdef MRCPP_HAS_MPI
    if (this->account_id >= 0) {
        MPI_Status status;
        int messages[message_size];
        messages[0] = GET_NEXTTASK;
        messages[1] = account_id;
        MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
        MPI_Recv(&nexttask, 1, MPI_INT, task_bank, 1, comm_bank, &status);
        return nexttask;
    }
#endif
    nexttask = this->task++;
    if (nexttask >= this->n_tasks) nexttask = -1;
    return nexttask;
}

void TaskManager::put_readytask(int i, int j) {
#ifdef MRCPP_HAS_MPI
    if (this->account_id < 0) return;
    MPI_Status status;
    int messages[message_size];
    messages[0] = PUT_READYTASK;
    messages[1] = account_id;
    messages[2] = i;
    messages[3] = j;
    MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
#endif
}

void TaskManager::del_readytask(int i, int j) {
#ifdef MRCPP_HAS_MPI
    if (this->account_id < 0) return;
    MPI_Status status;
    int messages[message_size];
    messages[0] = DEL_READYTASK;
    messages[1] = account_id;
    messages[2] = i;
    messages[3] = j;
    MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
#endif
}

std::vector<int> TaskManager::get_readytask(int i, int del) {
    std::vector<int> readytasks;
#ifdef MRCPP_HAS_MPI
    if (this->account_id < 0) return readytasks;
    MPI_Status status;
    int messages[message_size];
    messages[0] = GET_READYTASK;
    if (del == 1) messages[0] = GET_READYTASK_DEL;
    messages[1] = account_id;
    messages[2] = i;
    MPI_Send(messages, message_size, MPI_INT, task_bank, 0, comm_bank);
    int nready;
    MPI_Recv(&nready, 1, MPI_INT, task_bank, 844, comm_bank, &status);
    if (nready > 0) {
        readytasks.resize(nready);
        MPI_Recv(readytasks.data(), nready, MPI_INT, task_bank, 845, comm_bank, &status);
    }
#endif
    return readytasks;
}

} // namespace mrcpp
