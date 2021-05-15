/*
FastCDS -- Copyright (c) 2020, Xindi Zhang, Bohan Li, Shaowei Cai, et. al. dezhangxd@163.com

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <unordered_set>
#include <chrono>
#include <vector>
#include <algorithm>
#include <functional>
#include <stdint.h>


// some macros 
// used in the improved_1hop
// #define USE_RAND_IN_CONSTRUCTION
#define LOOKAHEAD_PICK_1HOP
// when fix_points, try not to use it because it is better! 
// but in our paper we use it because it is hard to experss in words.
#define USE_TIME_JUDGE
// #define ANALYZE_HDC
// #define COUNT_PICK_REMOVE
// using TABU/CC Method in NuCDS
// #define ORIGIN_CC_TABU_METHOD
// using pick/add method in NuCDS
// #define ORIGIN_PICK_REMOVE
// #define USE_SAFETY
// need to select which reasoning rules to use.
#define REASONING_MORE_PROB
// #define REASONING_FORCE
// #define RAND_BREAKING_TIE
#define USE_SAFETY_AGE

using namespace std;


// a heap
template<class T>
class My_heap;

// options
class Option;
class IntOption;
class DoubleOption;



class FastCDS{
    template<class T>
    friend class My_heap;
private:
    //control informations
    uint32_t        _opt_seed;
    uint8_t         _opt_verbosity;
    uint16_t        _opt_fix_type;
    uint32_t        _opt_bms_construct;
    uint32_t        _opt_bms_pick_add;
    uint32_t        _opt_bms_pick_remove;
    uint8_t         _opt_show_CDS;
    uint8_t         _opt_only_output_init_soln;
    uint32_t        _opt_reconstruct_gap;
    uint8_t         _opt_only_init;
    uint32_t        reconstruct_num;

    /* For evey restarts, use SUB\TBC iteratively.
     *  if step_this_turn < max_steps(1) || tmp_no_improve_steps < max_steps_noimpr(2), restart and switch.
     *  For (1) increase max_no_improve_steps, for(2) decrease.
     *  ‘base’ is the minimun unit(step) for change;
     *  if(2) best_cds_size < 5k, min_base_SUB*=10;
     *  if improve 10 steps & (1), max_steps_noimpr_SUB -= base;
     */
    uint32_t        restarts;
    uint64_t        steps;
    uint64_t        SUB_steps;
    uint32_t        steps_this_turn;
    uint32_t        no_improve_steps;
    uint32_t        max_steps_SUB;
    uint32_t        max_steps_TBC;
    uint32_t        tmp_max_steps_noimpr_SUB;
    uint32_t        tmp_max_steps_noimpr_TBC;
    uint32_t        min_base_SUB;
    uint32_t        min_base_TBC;
    uint32_t        try_step_SUB;
    uint32_t        try_step_TBC;
    bool            reconstruct_signal;
    uint8_t         ineffective_search_num;
    uint32_t        _opt_check_soln;
    uint32_t        _opt_only_SUB;
    uint32_t        _opt_smart_mv_multiedges;
    bool            TBC;


#ifdef RAND_BREAKING_TIE
    uint32_t        *break_tie_vec;
    uint32_t        break_tie_vec_sz;
#endif


    //basic informations
    uint32_t        v_nums;
    uint32_t        e_nums;
    uint32_t        max_degree;
    double          avg_degree;
    uint32_t        max_degree_vertex;
    uint32_t        *e_from,
                    *e_to;
    uint32_t        *degree;
    uint32_t        **adj;
    uint8_t         *is_grey;
    uint32_t        grey_vec_size;
    uint32_t        *grey_vec;
    uint32_t        *idx_in_grey_vec;
    uint32_t        white_vec_size;
    uint32_t        *white_vec;
    uint32_t        *idx_in_white_vec;
    uint32_t        undom_size;
    
#ifdef COUNT_PICK_REMOVE
    uint32_t        SUB_pick_num;
    uint32_t        TBC_pick_num;
    uint32_t        *ct_remove_SUB;
    uint32_t        *ct_remove_TBC;
#endif
    
    // tabu and CC
    uint64_t        *tabu_remove;
    uint64_t        *tabu_add;
    uint32_t        last_add;
    uint32_t        last_remove;
    uint8_t         *conf_change; //1:can rm/add, 0 fixed
    
    
    //Results
    uint32_t        best_cds_size;
    uint8_t         *best_cds;          //1: in cds ; 0: out of cds;
    uint32_t        cds_size;
    uint8_t         *cds;
    uint32_t        *cds_vec;
    uint32_t        *idx_in_cds_vec;
    double          best_cds_time;

    
    //for cut points && fix points
    uint32_t        *idx_visit;         // the visited neighbors index
    uint32_t        *dnf;               // the visit index in DFS.
    uint32_t        *low;               // the smallest dnf the node can reach.

    /* fix variables
    *   0 donot fix;
    *   1 fixed in S;
    *   2 fixed out S;
    */
    uint32_t        fixed_num;
    uint8_t         *fix_type;
    uint32_t        fix_1_vec_size;
    uint32_t        *fix_1_vec;
    uint8_t         *is_cut;
    uint32_t        cut_vec_size;
    uint32_t        *cut_vec;
    void            renew_cut();
    
    
    //used in SUB & TBC
    bool            increase_SUB_flag;
    uint32_t        root;
    uint32_t        *son_num;
    uint32_t        *leaf_vec;
    uint32_t        leaf_vec_size;
    uint32_t        *idx_in_leaf_vec;
    uint32_t        *father;
    uint8_t         *is_leaf;
    uint8_t         *seen;
    uint8_t         *seen_before_queue;
    void            construct_tree_random_based();
    void            construct_tree_random_based_for_big();
    void            construct_tree_hueristic_based();
    uint32_t        *best_father_cand;
    uint32_t        best_father_cand_size;
    
    // For subgraph in construct tree
    uint32_t        *subgraph_dom_num;
    uint32_t        *subgraph_dom_by;
    uint32_t        *subgraph_score;
    uint32_t        **subgraph_adj;
    uint8_t         *subgraph_is_grey;
    uint32_t        *subgraph_degree;
    My_heap<int32_t>*subgraph_heap;
    inline  void    set_subgraph_data_structure();
    inline  void    releases_subgraph_data_structure();
    void            construct_tree_subgraph_based();

    // some frequency information
    uint32_t        *pick_add_num;
    uint32_t        *pick_rm_num;
    uint32_t        *pick_swap_num;
    

    //different scores
    uint64_t        *time_stamp;
    /*  used for construction and cal the dom_num increasingly.
     *  grey <-> white:  white/grey nodes in N1[v] should +/- 1.
     *  black <-> grey： black nodes in N2[v] should +/- 1 if
     *      dom_num[]>1 related N1[v] . score[v] = -score[v]
     *
     */
    int32_t         *score;
    /* frequency, for each node u, freq[u] = 1;
     * for each step, node not in cds, freq[u]++;
     * for u in cds, freq_score = sum_freq(dom->undom)
     * for u not in cds, freq_score = sum_freq(undom->dom)
     */
    uint32_t        *freq;
    int32_t         *freq_score;
    /*  domination degree (dd) dd[u] = |N[G] in cds|
     *  safety_score[u] = -dd[u] for u in CDS
     *  safety_score[u] = dd[u] for u not in CDS
     *  is equal to +/-dom_num[u];
     *  !!!(do not use this score, use dom_num replace!!)
     */
    int32_t         *safety_score;
    // dominated by how many black nodes
    uint32_t        *dom_num;
    // dominated by which black nodes, if dom_num=1
    uint32_t        *dom_by;
    
#ifdef ANALYZE_HDC
    set<uint32_t>   before_SUB;
    set<uint32_t>   after_SUB;
    map<uint32_t,uint32_t>  TBC_add;
    map<uint32_t,uint32_t>  TBC_remove;
    map<uint32_t,uint32_t>  TBC_remove_last;
    map<uint32_t,uint32_t>  SUB_add;
    map<uint32_t,uint32_t>  SUB_remove;
#endif
    
    //choose node and update scores;
    uint8_t         _opt_for_small_ins;
    void            add_to_cds      (uint32_t vertex);
    void            remove_from_cds (uint32_t vertex);
    uint32_t        *pick_cand;
    uint32_t        pick_cand_size;
    uint32_t        pick_add();
    uint32_t        pick_add_small();
    uint32_t        pick_remove();
    uint32_t        pick_remove_small();
    void            qsort_by_score(uint32_t *a,uint32_t size);
    void            update_leaf_add(uint32_t vertex);
    void            update_leaf_remove(uint32_t vertex);
    
    
    //a stack
    uint32_t        my_stack_size;
    uint32_t        *my_stack;
    int32_t         *idx_in_my_stack;
    //a queue
    uint32_t        *my_queue;
    uint32_t        q_front;
    uint32_t        q_back;
    
    
    
    
    //runtime : do not take read&build time into account.
    double          last_time_point;
    double          _opt_cutoff_time;
    double          SUB_sum_time,TBC_sum_time;
    chrono::steady_clock::time_point start_time;
    double          getCpuTime();
    double          _opt_print_soln_gap;
    bool            should_print_soln;
    
    
    
    //memory
    void            alloc_memory();
    void            free_memory();
    
    
    //update_weight
    uint32_t        smooth_weight_ct;
    uint32_t        weight_threshold;
    uint32_t        total_weight;
    void            smooth_weight();
    
    
    //used in construction
    string          outer_cds_filename;
    void            reconstruct_with_upper_bound(uint32_t upper_bound);
    uint32_t        reconstruct_ct;
    
    //main functions for solving
    void            initialize();
    uint32_t        fix_points();
    //fix all vertices using a fast strategy.
    uint32_t        fix_points_fast();
    bool            check_soln();
    bool            check_tmp_soln_with_score();
    void            update_best_soln();
    void            update_freq_scores();
    
    
    
    
public:
    FastCDS();
    bool            build(string filename);
    uint32_t        construct_1hop();
    uint32_t        construct_1hop_bms();
    uint32_t        construct_by_outer();
    uint32_t        do_LS_SUB();
    uint32_t        do_LS_TBC();
    uint32_t        do_LS_TBC_without_HDC();
    bool            exec(string filename);
    bool            exec_only_TBC_plus();
    bool            exec_for_classical();
    bool            exec_for_massive();
    void            print_soln(bool need_check);
};



// some useful utils
// ---------------------------------------
template<class T>
class My_heap{
public:
    //my heap
    bool            cmp_using_function;
    function<bool(uint32_t,uint32_t)>     cmp;
    uint32_t        mem_size;
    T const         *heap_score;
    uint32_t        my_heap_size;
    uint32_t        *my_heap;
    int32_t         *idx_in_my_heap;
    inline int32_t  father      (int32_t idx) const;
    inline int32_t  left_son    (int32_t idx) const;
    inline int32_t  right_son   (int32_t idx) const;
    inline void     conduct_up  (int32_t idx);
    inline void     conduct_down(int32_t idx);
    inline void     remove_node (uint32_t vertex);
    inline void     add_node    (uint32_t vertex);
    inline uint32_t pop_top     ();
    My_heap(uint32_t mem_size, const T *heap_score);
    My_heap(uint32_t mem_size, function<bool(uint32_t,uint32_t)>fun);
    ~My_heap();
};




class Option{
protected:
    const char* name;
    const char* help;
    Option(const char* name,
           const char* help):
    name(name),help(help){
        get_options().push_back(this);}
    static vector<Option*>& get_options(){
        static vector<Option*> options;
        return options;}
public:
    const char*   get_name(){return name;}
    const char*   get_help(){return help;}
    virtual string get_default()=0;
    virtual string get_value()=0;
    virtual void  parse(string value)=0;
    
    friend  void parse_options(int argc, char *argv[]);
    friend  void print_uses_and_exit();
    friend  void print_parameters();
};


class IntOption : public Option{
protected:
    int32_t value,default_value;
public:
    IntOption(const char* name, const char* help, int32_t default_value):
    Option(name,help),value(default_value){value = default_value;}
    virtual void parse(string s_value){
        stringstream ss(s_value);ss>>value;}
    virtual string get_default(){
        stringstream ss;ss<<default_value;return ss.str();}
    virtual string get_value(){
        stringstream ss;ss<<value;return ss.str();}
    operator      int32_t   (void) const {return value; }
    operator      int32_t&  (void)       {return value; }
    IntOption&    operator= (int32_t x)  {value = x; return *this;}
};
class DoubleOption : public Option{
protected:
    double value,default_value;
public:
    DoubleOption(const char* name, const char* help, double default_value):
    Option(name,help),value(default_value){value = default_value;}
    virtual void parse(string s_value){
        stringstream ss(s_value);ss>>value;}
    virtual string get_default(){
        stringstream ss;ss<<default_value;return ss.str();}
    virtual string get_value(){
        stringstream ss;ss<<value;return ss.str();}
    operator      double   (void) const {return value; }
    operator      double&  (void)       {return value; }
    DoubleOption&    operator= (double x)  {value = x; return *this;}
};
//class StringOption : public Option{
//protected:
//    const char* value;
//    const char* default_value;
//public:
//    StringOption(const char* name, const char* help, const char* default_value=NULL):
//    Option(name,help),value(default_value){value = default_value;}
//    virtual void parse(string s_value){value=s_value.c_str();cout<<value<<endl;}
//    virtual string get_default(){
//        if(default_value==NULL) return string("null");else return string(default_value);}
//    virtual string get_value(){
//        if(value==NULL) return string("null");else return string(value);}
//    operator const char*   (void) const {return value;cout<<value<<endl; }
//        operator const char*&  (void)       {return value;cout<<value<<"2"<<endl; }
//    StringOption&    operator= (const char* x)  {value = x;cout<<value<<"3"<<endl; return *this;}
//};
class StringOption : public Option{
protected:
    string value;
    string default_value;
public:
    StringOption(const char* name, const char* help, const char* default_value=""):
    Option(name,help),value(default_value){value = default_value;}
    virtual void parse(string s_value){value=s_value.c_str();}
    virtual string get_default(){return default_value;}
    virtual string get_value(){return value;}
    operator string   (void) const {return value;cout<<value<<endl; }
    operator string&  (void)       {return value;}
    StringOption&    operator= (const char* x)  {value = x; return *this;}
};

void print_uses_and_exit(){
    for(auto& option : Option::get_options()){
        cout<<"--"<<option->get_name()<<" (";
        cout<<"default:\t"<<option->get_default()<<")"<<endl;
        cout<<option->get_help()<<endl;
    }
}
void print_parameters(){
    for(auto& option : Option::get_options()){
        cout<<"c --PARM--" <<option->get_name()<<"("<<option->get_value()<<")"<<endl;
    }
}
void parse_options(int argc, char** argv){
    int filename_ct = 0;
    for(int i=1;i<argc;++i){
        bool find_same = false;
        string s_para(argv[i]);
        if(s_para==string("help") || s_para==string("-h") || s_para==string("--help")){
            cout<<"c use ./fastCDS [--option_name=value] <input-ins>"<<endl;
            print_uses_and_exit();
        }
        if(! (s_para.size()>2 && s_para[0]=='-' && s_para[1]=='-') ){
            filename_ct++;
            argv[1] = argv[i];
            continue;
        }

        if(s_para.size()>4 && filename_ct<2 && s_para.find('=')>1){
            string subs = s_para.substr(2,s_para.find('=')-2);
            for(auto& option : Option::get_options()){
                if(string(option->get_name()) == subs){
                    string value_s=s_para.substr(s_para.find('=')+1);
                    option->parse(value_s);
                    find_same = true;
                    break;
                }
            }
        }
        if(!find_same){
            cout<<"c --PARSE ERROR--"<<endl;
            cout<<"c use ./fastCDS [--option_name=value] <input-ins>"<<endl;
            exit(404);
        }
    }
    if(filename_ct!=1){
        cout<<"c --PARSE ERROR--"<<endl;
        cout<<"c use ./fastCDS [--option_name=value] <input-ins>"<<endl;
        exit(404);
    }
}

