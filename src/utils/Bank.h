#pragma once

#include "ComplexFunction.h"
#include "parallel.h"
#include "trees/NodeIndex.h"

namespace mrcpp {

using namespace mpi;

struct deposit {
    ComplexFunction *orb;
    double *data; // for pure data arrays
    bool hasdata;
    int datasize;
    int id = -1; // to identify what is deposited
    int source;  // mpi rank from the source of the data
};

struct queue_struct {
    int id;
    std::vector<int> clients;
};

enum {
    // (the values are used to interpret error messages)
    CLOSE_BANK,              // 0
    CLEAR_BANK,              // 1
    NEW_ACCOUNT,             // 2
    CLOSE_ACCOUNT,           // 3
    GET_ORBITAL,             // 4
    GET_FUNCTION_AND_WAIT,   // 5
    GET_FUNCTION_AND_DELETE, // 6
    SAVE_ORBITAL,            // 7
    GET_FUNCTION,            // 8
    SAVE_FUNCTION,           // 9
    GET_DATA,                // 10
    SAVE_DATA,               // 11
    SAVE_NODEDATA,           // 12
    GET_NODEDATA,            // 13
    GET_NODEBLOCK,           // 14
    GET_ORBBLOCK,            // 15
    CLEAR_BLOCKS,            // 16
    GET_MAXTOTDATA,          // 17
    GET_TOTDATA,             // 18
    INIT_TASKS,              // 19
    GET_NEXTTASK,            // 20
    PUT_READYTASK,           // 21
    DEL_READYTASK,           // 22
    GET_READYTASK,           // 23
    GET_READYTASK_DEL,       // 24
};

class Bank {
public:
    Bank() = default;
    ~Bank();
    void open();
    void close();
    int get_maxtotalsize();
    std::vector<int> get_totalsize();

private:
    friend class BankAccount;
    friend class TaskManager;

    // used by BankAccount
    int openAccount(int iclient, MPI_Comm comm);
    int clearAccount(int account, int iclient, MPI_Comm comm); // closes and open fresh account
    void closeAccount(int account_id);                         // remove the account

    // used by TaskManager;
    int openTaskManager(int ntasks, int iclient, MPI_Comm comm);
    void closeTaskManager(int account_id);

    // used internally by Bank;
    void clear_bank();
    void remove_account(int account); // remove the content and the account

    long long totcurrentsize = 0ll;                     // number of kB used by all accounts
    std::vector<int> accounts;                          // open bank accounts
    std::map<int, std::vector<deposit> *> get_deposits; // gives deposits of an account
    std::map<int, std::map<int, int> *> get_id2ix;
    std::map<int, std::map<int, int> *> get_id2qu;
    std::map<int, std::vector<queue_struct> *> get_queue;            // gives deposits of an account
    std::map<int, std::map<int, std::vector<int>> *> get_readytasks; // used by task manager
    std::map<int, long long> currentsize;                            // total deposited data size (without containers)
    long long maxsize = 0;                                           // max total deposited data size (without containers)
};

class BankAccount {
public:
    BankAccount(int iclient = wrk_rank, MPI_Comm comm = comm_wrk);
    ~BankAccount();
    int account_id = -1;
    void clear(int i = wrk_rank, MPI_Comm comm = comm_wrk);
    //    int put_orb(int id, ComplexFunction &orb);
    //    int get_orb(int id, ComplexFunction &orb, int wait = 0);
    int get_func_del(int id, ComplexFunction &orb);
    int put_func(int id, ComplexFunction &func);
    int get_func(int id, ComplexFunction &func, int wait = 0);
    int put_data(int id, int size, double *data);
    int get_data(int id, int size, double *data);
    int put_data(NodeIndex<3> nIdx, int size, double *data);
    int get_data(NodeIndex<3> nIdx, int size, double *data);
    int put_nodedata(int id, int nodeid, int size, double *data);
    int get_nodedata(int id, int nodeid, int size, double *data, std::vector<int> &idVec);
    int get_nodeblock(int nodeid, double *data, std::vector<int> &idVec);
    int get_orbblock(int orbid, double *&data, std::vector<int> &nodeidVec, int bankstart);
};

class TaskManager {
public:
    TaskManager(int ntasks, int iclient = wrk_rank, MPI_Comm comm = comm_wrk);
    ~TaskManager();
    int next_task();
    void put_readytask(int i, int j);
    void del_readytask(int i, int j);
    std::vector<int> get_readytask(int i, int del);
    int account_id = -1;
    int task = 0;    // used in serial case only
    int n_tasks = 0; // used in serial case only
};

int const message_size = 7;

} // namespace mrcpp
