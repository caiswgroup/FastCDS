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



#include "fastCDS.hpp"

// options
static IntOption opt_fix_type
 ("fix_type","\
    0       donnot use fix method;\n\
    1b      fixed cutPoint;\n\
    10b     fixed 0 with 2,3-clique;\n\
    100b    fixed all k-clique;\n\
    1000b   adjust atomatically;\n\
    (first 12 bits means the max k if 1000b is the last 4 bits);",0x0008
  );
static IntOption opt_bms_construct
 ("bms_construct","\
    0       use heap to construct;\n\
    >0      BMS pick num;",0
  );
static IntOption opt_verbosity
 ("verbosity","\
    0       print soln only;\n\
    1       print baisc info;\n\
    2       all",1
  );
static IntOption opt_seed("seed","",0);
static DoubleOption opt_cutoff_time("cutoff_time","",3600);
static DoubleOption opt_print_soln_gap("print_gap","",1000);
static IntOption opt_bms_pick_add("bms_pick_add","",10);
static IntOption opt_bms_pick_remove("bms_pick_remove","",100);
static IntOption opt_show_CDS("show_CDS","bool",0);
static StringOption opt_build_outer("build_by","","null");
static IntOption opt_check_soln("check_soln","bool",0);
static IntOption opt_only_SUB("only_SUB","bool",0);
static IntOption opt_smart_mv_multiedges("smart_mv_multiedges","bool",1);
static IntOption opt_only_output_init_soln("only_output_init_soln","bool",0);
static IntOption opt_for_small_ins("for_small_ins","bool",0);
static IntOption opt_reconstruct_gap("reconstruct_tree","",500);
static IntOption opt_only_init("only_init","bool",0);

FastCDS::FastCDS():
 _opt_verbosity(opt_verbosity)
,_opt_fix_type(opt_fix_type)
,_opt_bms_construct(opt_bms_construct)
,_opt_seed(opt_seed)
,_opt_cutoff_time(opt_cutoff_time)
,_opt_print_soln_gap(opt_print_soln_gap)
,_opt_bms_pick_add(opt_bms_pick_add)
,_opt_bms_pick_remove(opt_bms_pick_remove)
,outer_cds_filename(opt_build_outer.get_value())
,_opt_show_CDS(opt_show_CDS)
,_opt_check_soln(opt_check_soln)
,_opt_only_SUB(opt_only_SUB)
,_opt_smart_mv_multiedges(opt_smart_mv_multiedges)
,_opt_only_output_init_soln(opt_only_output_init_soln)
,_opt_for_small_ins(opt_for_small_ins)
,_opt_reconstruct_gap(opt_reconstruct_gap)
,_opt_only_init(opt_only_init)
{
    best_cds_time               = (double)0;
    best_cds_size               = INT32_MAX;
    fixed_num                   = 0;
    steps                       = 0;
    restarts                    = 0;
    min_base_SUB                = 100;
    min_base_TBC                = 100*1000;
    max_steps_SUB               = 10 * min_base_SUB;
    max_steps_TBC               = 10 * min_base_TBC;
    tmp_max_steps_noimpr_SUB    = 2  * min_base_SUB;
    tmp_max_steps_noimpr_TBC    = 2  * min_base_TBC;
    try_step_SUB                = 10;
    try_step_TBC                = 10000;
// #ifndef ANALYZE_HDC
    if(_opt_for_small_ins==0){
        min_base_SUB            = 100;
        max_steps_SUB           = 2  * min_base_SUB;
        tmp_max_steps_noimpr_SUB= 1  * min_base_SUB;
        try_step_SUB            = 10;
    }
// #endif
    last_time_point             = 0;
    increase_SUB_flag           = true;
    fix_1_vec_size              = 0;
    undom_size                  = 0; //after initailized.
    cut_vec_size                = 0;
    should_print_soln           = false;
    leaf_vec_size               = 0;
    SUB_steps                   = 0;
    smooth_weight_ct            = 0;
    reconstruct_ct              = 0;
    last_add                    = 0;
    last_remove                 = 0;
    reconstruct_num             = 0;

    srand(_opt_seed);
    if(_opt_verbosity>1) print_parameters();
}


void FastCDS::alloc_memory(){
    uint32_t v_mems     = v_nums + 1;
    uint32_t e_mems     = e_nums + 1;
    e_from          =   new uint32_t[e_mems];
    e_to            =   new uint32_t[e_mems];
    degree          =   new uint32_t[v_mems];
    adj             =   new uint32_t*[v_mems];
    idx_visit       =   new uint32_t[v_mems];
    dnf             =   new uint32_t[v_mems];
    low             =   new uint32_t[v_mems];
    fix_type        =   new uint8_t[v_mems];
    my_stack        =   new uint32_t[v_mems];
    time_stamp      =   new uint64_t[v_mems];
    score           =   new int32_t[v_mems];
    freq            =   new uint32_t[v_mems];
    freq_score      =   new int32_t[v_mems];
    safety_score    =   new int32_t[v_mems];
    best_cds        =   new uint8_t[v_mems];
    cds             =   new uint8_t[v_mems];
    is_grey         =   new uint8_t[v_mems];
    idx_in_my_stack =   new int32_t[v_mems];
    fix_1_vec       =   new uint32_t[v_mems];
    cds_vec         =   new uint32_t[v_mems];
    idx_in_cds_vec  =   new uint32_t[v_mems];
    dom_num         =   new uint32_t[v_mems];
    dom_by          =   new uint32_t[v_mems];
    cut_vec         =   new uint32_t[v_mems];
    is_cut          =   new uint8_t[v_mems];
    tabu_remove     =   new uint64_t[v_mems];
    tabu_add        =   new uint64_t[v_mems];
    conf_change     =   new uint8_t[v_mems];
    grey_vec        =   new uint32_t[v_mems];
    idx_in_grey_vec =   new uint32_t[v_mems];
    white_vec       =   new uint32_t[v_mems];
    idx_in_white_vec=   new uint32_t[v_mems];
    leaf_vec        =   new uint32_t[v_mems];
    is_leaf         =   new uint8_t[v_mems];
    idx_in_leaf_vec =   new uint32_t[v_mems];
    father          =   new uint32_t[v_mems];
    son_num         =   new uint32_t[v_mems];
    my_queue        =   new uint32_t[v_mems+1];
    seen            =   new uint8_t[v_mems];
    pick_add_num    =   new uint32_t[v_mems];
    pick_rm_num     =   new uint32_t[v_mems];
    pick_swap_num   =   new uint32_t[v_mems];
    best_father_cand=   new uint32_t[v_mems];
    seen_before_queue=  new uint8_t[v_mems];

#ifdef RAND_BREAKING_TIE
    break_tie_vec   =   new uint32_t[v_mems];
#endif
#ifdef COUNT_PICK_REMOVE
    ct_remove_TBC   = new uint32_t[v_mems];
    ct_remove_SUB   = new uint32_t[v_mems];
    fill_n(ct_remove_TBC,v_mems,0);
    fill_n(ct_remove_SUB,v_mems,0);
    SUB_pick_num = TBC_pick_num = 0;
#endif

    fill_n(seen_before_queue,v_mems, 0);
    fill_n(is_leaf,     v_mems, 0);
    fill_n(tabu_remove, v_mems, 0);
    fill_n(tabu_add, v_mems, 0);
    fill_n(conf_change, v_mems, 1);
    fill_n(dom_num,     v_mems, 0);
    fill_n(freq,        v_mems, 1);
    fill_n(degree,      v_mems, 0);
    fill_n(time_stamp,  v_mems, 0);
    fill_n(fix_type,    v_mems, 0);
    fill_n(is_cut,      v_mems, 0);
    fill_n(pick_add_num,v_mems, 0);
    fill_n(pick_rm_num, v_mems, 0);
    fill_n(pick_swap_num,v_mems, 0);
    // initalize some variables
    if(_opt_for_small_ins==1){
        weight_threshold    = (v_nums+1)*100;
    }else{
        weight_threshold    = (v_nums+1)*1000;
    }
    total_weight        = v_nums;
    time_stamp[0]       = INT64_MAX;
}

void FastCDS::free_memory(){
    delete [] e_from;
    delete [] e_to;
    delete [] degree;
    delete [] fix_type;
    delete [] idx_visit;
    delete [] time_stamp;
    delete [] low;
    delete [] dnf;
    delete [] my_stack;
    delete [] score;
    delete [] is_grey;
    delete [] idx_in_my_stack;
    delete [] freq;
    delete [] freq_score;
    delete [] safety_score;
    delete [] fix_1_vec;
    delete [] cds_vec;
    delete [] idx_in_cds_vec;
    delete [] dom_num;
    delete [] dom_by;
    delete [] is_cut;
    delete [] cut_vec;
    delete [] pick_cand;
    delete [] tabu_remove;
    delete [] conf_change;
    delete [] grey_vec;
    delete [] idx_in_grey_vec;
    delete [] white_vec;
    delete [] idx_in_white_vec;
    delete [] leaf_vec;
    delete [] is_leaf;
    delete [] father;
    delete [] idx_in_leaf_vec;
    delete [] son_num;
    delete [] my_queue;
    delete [] seen;
    delete [] pick_add_num;
    delete [] pick_rm_num;
    delete [] pick_swap_num;
    delete [] best_father_cand;
    delete [] seen_before_queue;
    for(uint32_t i=1;i<=v_nums;++i) delete[] adj[i];
    delete [] adj;
    if(_opt_for_small_ins!=1){
        releases_subgraph_data_structure();
    }
#ifdef RAND_BREAKING_TIE
    delete [] break_tie_vec;
#endif
#ifdef COUNT_PICK_REMOVE
    delete [] ct_remove_TBC;
    delete [] ct_remove_SUB;
#endif
}




//bsic fucntions
//===============================
double FastCDS::getCpuTime(){
    chrono::steady_clock::time_point now_time = chrono::steady_clock::now();;
    chrono::duration<double> duration = now_time - start_time;
    return duration.count();
}

bool FastCDS::check_soln(){
    //check cds size
    uint32_t size = 0;
    for(uint32_t idx=0;idx<=v_nums;++idx){
        if(best_cds[idx] == 1) size++;
    }
    if(size!=best_cds_size) return false;
    
    //check domination
    uint8_t *tmp_dom      = new uint8_t[v_nums+1];
    fill_n(tmp_dom,v_nums+1,0);
    for(uint32_t idx=1;idx<=v_nums;++idx){
        if(best_cds[idx] == 1){
            tmp_dom[idx] = 1;
            for(uint32_t n1_idx=0;n1_idx<degree[idx];++n1_idx){
                tmp_dom[adj[idx][n1_idx]] = 1;
            }
        }
    }
    bool result = true;
    for(uint32_t idx=1;idx<=v_nums;++idx){
        if(tmp_dom[idx]==0){
            result = false;
            break;
        }
    }
    delete [] tmp_dom;
    if(result == false) return result;
    
    
    //check connection
    uint8_t *tmp_soln = new uint8_t[v_nums+1];
    copy(best_cds,best_cds+v_nums+1,tmp_soln);
    my_stack_size = 0;
    for(uint32_t idx=1;idx<v_nums;++idx){
        if(tmp_soln[idx] == 1){
            my_stack[my_stack_size++] = idx;
            tmp_soln[idx] = 0;
            break;
        }
    }
    while(my_stack_size>0){
        uint32_t top = my_stack[--my_stack_size];
        for(uint32_t n1_idx=0;n1_idx<degree[top];++n1_idx){
            uint32_t n1 = adj[top][n1_idx];
            if(tmp_soln[n1] == 1){
                my_stack[my_stack_size++] = n1;
                tmp_soln[n1] = 0;
            }
        }
    }
    result = true;
    for(uint32_t idx=1;idx<=v_nums;++idx){
        if(tmp_soln[idx]==1){
            result = false;
            break;
        }
    }
    delete [] tmp_soln;
    return result;
}


bool FastCDS::check_tmp_soln_with_score(){
    // check cds\white\grey idx & size;
    uint32_t size_cds = 0;
    for(uint32_t node=1;node<=v_nums;++node){if(cds[node]==1) size_cds++;}
    if(size_cds!=cds_size) cout<<"c cds size wrong right="<<size_cds<<" wrong="<<cds_size<<endl;
    for(uint32_t idx=0;idx<cds_size;++idx){
        if(idx_in_cds_vec[cds_vec[idx]] != idx){cout<<"c cds idx wrong"<<endl;}
        if(cds[cds_vec[idx]]==0){cout<<"c not in cds"<<endl;}
    }
    uint32_t size_grey = 0;
    for(uint32_t node=1;node<=v_nums;++node){if(is_grey[node]==1) size_grey++;}
    if(size_grey!=grey_vec_size) cout<<"c grey size wrong"<<endl;
    for(uint32_t idx=0;idx<grey_vec_size;++idx){
        if(idx_in_grey_vec[grey_vec[idx]] != idx){cout<<"c grey idx wrong"<<endl;}
        if(is_grey[grey_vec[idx]]==0){cout<<"c not in grey"<<endl;}
    }
    uint32_t size_white = 0;
    for(uint32_t node=1;node<=v_nums;++node){if(cds[node]==0&&is_grey[node]==0) size_white++;}
    if(size_white!=white_vec_size) cout<<"c white size wrong"<<endl;
    for(uint32_t idx=0;idx<white_vec_size;++idx){
        if(idx_in_white_vec[white_vec[idx]] != idx){cout<<"c white idx wrong"<<endl;}
        if(cds[white_vec[idx]]==1 || is_grey[white_vec[idx]]==1){cout<<"c not in white"<<endl;}
    }
    if(size_cds+size_grey+size_white!=v_nums){cout<<"c sum size error"<<endl;}
    if(undom_size != white_vec_size){cout<<"c undom size error"<<endl;}
    //check score
    for(uint32_t node=1;node<=v_nums;++node){
        int32_t tmp_score=0;
        int32_t tmp_freq_score=0;
        uint32_t tmp_dom_num=0;
        uint32_t tmp_dom_by=0;
        if(cds[node]==1){
            tmp_dom_num++;
            tmp_dom_by = node;
            if(dom_num[node]==1){
                tmp_score--;
                tmp_freq_score-=freq[node];
            }
            for(uint32_t n1_idx=0;n1_idx<degree[node];++n1_idx){
                uint32_t n1 = adj[node][n1_idx];
                if(dom_num[n1]==1){
                    tmp_score--;
                    tmp_freq_score-=freq[n1];
                }
            }
        }else{
            if(dom_num[node]==0){
                tmp_freq_score += freq[node];
                tmp_score += 1;
            }
            for(uint32_t n1_idx=0;n1_idx<degree[node];++n1_idx){
                uint32_t n1 = adj[node][n1_idx];
                if(dom_num[n1]==0){
                    tmp_freq_score += freq[n1];
                    tmp_score += 1;
                }
            }
        }
        if(tmp_score != score[node]){cout<<"c "<<node<<" score wrong ,right="<<tmp_score<<" score"<<score[node]<<endl;}
        if(tmp_freq_score != freq_score[node]){cout<<"c "<<node<<" freq score wrong ,right="<<tmp_freq_score<<" fscore"<<freq_score[node]<<endl;}
        
        for(uint32_t n1_idx=0;n1_idx<degree[node];++n1_idx){
            uint32_t n1 = adj[node][n1_idx];
            if(cds[n1]==1){
                tmp_dom_by = n1;
                tmp_dom_num++;
            }
        }
        if(tmp_dom_num != dom_num[node]){
            cout<<"c dom num wrong"<<endl;
        }else if(tmp_dom_by!=dom_by[node] && tmp_dom_num==1){
            cout<<"c dom by wrong:"<<node<<" by "<<tmp_dom_by<<endl;
        }
    }
    //check cds connectivity;
//    for(uint32_t idx=1;idx<=v_nums;++idx){
//        if(cds[idx]==1) cout<<idx<<" ";
//    }cout<<endl;
    uint32_t *tmp_stack = new uint32_t[v_nums+1];
    uint8_t *seen = new uint8_t[v_nums+1];
    fill_n(seen,v_nums+1,0);
    uint32_t  tmp_stack_size = 1;
    tmp_stack[0] = cds_vec[0];
    seen[tmp_stack[0]]=1;
    while(tmp_stack_size>0){
        uint32_t top = tmp_stack[--tmp_stack_size];
        for(uint32_t idx=0;idx<degree[top];++idx){
            uint32_t n1 = adj[top][idx];
            if(cds[n1]==1 && seen[n1]==0){
                seen[n1] = 1;
                tmp_stack[tmp_stack_size++] = n1;
            }
        }
    }
    for(uint32_t idx=1;idx<=v_nums;++idx){
        if(cds[idx]==1 && seen[idx]==0){cout<<"c not connected"<<endl;}
    }
    delete [] tmp_stack;
    delete [] seen;
    
    //check the correctness of color assign
    for(uint32_t node=1;node<=v_nums;++node){
        if(cds[node]==1){}
        else{
            if(dom_num[node]>0){
                if(is_grey[node]==0){cout<<"c grey node assign error"<<endl;}
            }else{
                if(is_grey[node]==1){cout<<"c white node assign error"<<endl;}
            }
        }
    }
//    cout<<"white:"<<white_vec_size<<"\t";for(uint32_t idx=0;idx<white_vec_size;++idx){cout<<white_vec[idx]<<" ";}cout<<endl;
//    cout<<"grey:"<<grey_vec_size<<"\t";for(uint32_t idx=0;idx<grey_vec_size;++idx){cout<<grey_vec[idx]<<" ";}cout<<endl;
//    cout<<"cds:"<<cds_size<<"\t";for(uint32_t idx=0;idx<cds_size;++idx){cout<<cds_vec[idx]<<" ";}cout<<endl;

    
    return true;
}



void FastCDS::print_soln(bool need_check){
    if(_opt_show_CDS>0 && need_check && !check_soln()){
        cout<<"c --BUG-- The soln result is wrong!"<<endl;
    }
    cout<<"s "<<best_cds_size<<" "<<best_cds_time<<endl;
    if(_opt_show_CDS > 0){
        cout<<"v ";
        for(uint32_t i=1;i<=v_nums;++i){
            if(best_cds[i] == 1) cout<<i<<" ";
        }cout<<endl;
    }
}


void FastCDS::update_best_soln(){
    if(cds_size < best_cds_size ){
        double tmp_time = getCpuTime();
        if(tmp_time > _opt_cutoff_time) return;
        best_cds_size   = cds_size;
        best_cds_time   = tmp_time;
        if(_opt_show_CDS>0) copy(cds,cds+v_nums+1,best_cds);
        if(should_print_soln){cout<<"o "<<best_cds_size<<" "<<best_cds_time<<" "<<steps<<" "<<restarts<<endl;}
        if(_opt_check_soln!=0){check_tmp_soln_with_score();}
    }
}


//heap related
//-------------------------------
template<class T>
My_heap<T>::My_heap(uint32_t mem_size, const T *heap_score):
mem_size(mem_size),
heap_score(heap_score)
{   
    cmp_using_function  = false;
    my_heap_size        = 0;
    uint32_t v_mems     = mem_size+1;
    my_heap         =   new uint32_t[v_mems];
    idx_in_my_heap  =   new int32_t[v_mems];
    fill_n(idx_in_my_heap,v_mems,-1);
}
template<class T>
My_heap<T>::My_heap(uint32_t mem_size, function<bool(uint32_t,uint32_t)> fun):
mem_size(mem_size),
cmp(fun)
{
    cmp_using_function  = true;
    my_heap_size        = 0;
    uint32_t v_mems     = mem_size+1;
    my_heap         =   new uint32_t[v_mems];
    idx_in_my_heap  =   new int32_t[v_mems];
    fill_n(idx_in_my_heap,v_mems,-1);
}
// template<class T>
// My_heap<T>::My_heap(uint32_t mem_size, int use_function_symbol):
// mem_size(mem_size)
// {
//     cmp_using_function  = true;
//     my_heap_size        = 0;
//     uint32_t v_mems     = mem_size+1;
//     my_heap         =   new uint32_t[v_mems];
//     idx_in_my_heap  =   new int32_t[v_mems];
//     fill_n(idx_in_my_heap,v_mems,-1);
// }
template<class T>
My_heap<T>::~My_heap(){
    delete [] my_heap;
    delete [] idx_in_my_heap;
}
template<class T>
int32_t My_heap<T>::father      (int32_t idx) const {return (idx-1)/2;}
template<class T>
int32_t My_heap<T>::left_son    (int32_t idx) const {return idx*2+1;}
template<class T>
int32_t My_heap<T>::right_son   (int32_t idx) const {return idx*2+2;}
template<class T>
void My_heap<T>::conduct_up  (int32_t idx){
    if(cmp_using_function){
        while(idx>0){
            int32_t f          = father(idx);
            uint32_t v_idx      = my_heap[idx];
            uint32_t v_father   = my_heap[f];
            if(cmp(v_father,v_idx)){
                swap(my_heap[f],my_heap[idx]);
                swap(idx_in_my_heap[v_father],idx_in_my_heap[v_idx]);
                idx = f;
            }else return;
        }
    }else{
        while(idx>0){
            int32_t f          = father(idx);
            uint32_t v_idx      = my_heap[idx];
            uint32_t v_father   = my_heap[f];
            if(heap_score[v_idx]>heap_score[v_father]){
                swap(my_heap[f],my_heap[idx]);
                swap(idx_in_my_heap[v_father],idx_in_my_heap[v_idx]);
                idx = f;
            }else return;
        }
    }
}
template<class T>
void My_heap<T>::conduct_down(int32_t idx){
    int32_t left = left_son(idx);
    int32_t right= right_son(idx);
    if(cmp_using_function){
        while(left < my_heap_size){
            int32_t maxson   = left;
            if(right < my_heap_size && cmp(my_heap[left],my_heap[right]))
                maxson = right;
            int32_t v_maxson = my_heap[maxson];
            uint32_t v_idx    = my_heap[idx];
            if(cmp(v_idx,v_maxson)){
                swap(my_heap[maxson],my_heap[idx]);
                swap(idx_in_my_heap[v_maxson],idx_in_my_heap[v_idx]);
            }else return;
            idx     = maxson;
            left    = left_son(idx);
            right   = right_son(idx);
        }
    }else{
        while(left < my_heap_size){
            int32_t maxson   = left;
            if(right < my_heap_size && heap_score[my_heap[left]] < heap_score[my_heap[right]])
                maxson = right;
            int32_t v_maxson = my_heap[maxson];
            uint32_t v_idx    = my_heap[idx];
            if(heap_score[v_maxson] > heap_score[v_idx]){
                swap(my_heap[maxson],my_heap[idx]);
                swap(idx_in_my_heap[v_maxson],idx_in_my_heap[v_idx]);
            }else return;
            idx     = maxson;
            left    = left_son(idx);
            right   = right_son(idx);
        }
    }
}
template<class T>
void My_heap<T>::remove_node (uint32_t vertex){
    int32_t idx_vertex = idx_in_my_heap[vertex];
    uint32_t last_vertex = my_heap[--my_heap_size];
    idx_in_my_heap[last_vertex] = idx_vertex;
    my_heap[idx_vertex] = last_vertex;
    idx_in_my_heap[vertex] = -1;
    if(cmp_using_function){
        if(cmp(vertex,last_vertex))
            conduct_up(idx_vertex);
        else
            conduct_down(idx_vertex);
    }else{
        if(heap_score[last_vertex] > heap_score[vertex])
            conduct_up(idx_vertex);
        else
            conduct_down(idx_vertex);
    }
}
template<class T>
void My_heap<T>::add_node    (uint32_t vertex){
    idx_in_my_heap[vertex] = my_heap_size;
    my_heap[my_heap_size] = vertex;
    conduct_up(my_heap_size++);
}
template<class T>
uint32_t My_heap<T>::pop_top(){
    uint32_t top = my_heap[0];
    uint32_t last = my_heap[--my_heap_size];
    idx_in_my_heap[top] = -1;
    idx_in_my_heap[last] = 0;
    my_heap[0] = last;
    conduct_down(0);
    return top;
}
// template<class T>
// void My_heap<T>::set_callback_func(function<bool(uint32_t,uint32_t)>fun){
//     cmp = fun;
// }

// For initialize
//========================================================
bool FastCDS::build(string filename){
    ifstream fin(filename);
    if(!fin)return false;
    string line, tmp;
    stringstream buffer;
    buffer << fin.rdbuf();
    fin.close();
    start_time = chrono::steady_clock::now();
    stringstream ss;
    getline(buffer,line);
    stringstream ss_info(line);
    if(line[0] == 'p')
        ss_info>>tmp>>tmp>>v_nums>>e_nums;
    else
        ss_info>>v_nums>>e_nums;

    alloc_memory();
    
    for(uint32_t i=0;i<e_nums;++i){
        getline(buffer,line);
        stringstream ss(line);
        if(line[0]=='p')
            ss>>tmp>>e_from[i]>>e_to[i];
        else
            ss>>e_from[i]>>e_to[i];
        ++degree[e_from[i]];
        ++degree[e_to[i]];
    }
    max_degree = 0; max_degree_vertex = -1;
    for(uint32_t i=1;i<=v_nums;++i) {
        adj[i] = new uint32_t[degree[i]+1];
    }
    uint32_t  *degree_tmp = new uint32_t[v_nums+2];
    fill_n(degree_tmp,v_nums+2,0);
    uint32_t e_nums_origin = e_nums;
    for(uint32_t i=0;i<e_nums_origin;++i){
        uint32_t u = e_from[i];
        uint32_t v = e_to[i];
        if(u == v) {e_nums--;continue;}
        bool dumplicate = false;
        if(_opt_smart_mv_multiedges==0 || (_opt_smart_mv_multiedges!=0 && v_nums<1000) ){
            for(uint32_t adj_idx=0;adj_idx<degree_tmp[u];++adj_idx){
                if(adj[u][adj_idx] == v) {
                    e_nums--;
                    dumplicate=true;
                    break;
                }
            }
        }
        if(dumplicate) continue;
        adj[u][degree_tmp[u]++] = v;
        adj[v][degree_tmp[v]++] = u;
    }
    for(uint32_t i=1;i<=v_nums;++i) {
        if(degree[i] > max_degree){
            max_degree = degree[i];
            max_degree_vertex = i;
        }
    }

    copy(degree_tmp, degree_tmp+v_nums+1, degree);
    delete [] degree_tmp;
    avg_degree = (double)2*e_nums/v_nums;

    uint32_t pick_cand_size;
    if(_opt_for_small_ins == 1){
        pick_cand_size = max_degree*v_nums;
    }else{
        pick_cand_size = v_nums + 1;
        pick_cand_size = max(pick_cand_size,max_degree*_opt_bms_pick_add);
        pick_cand_size = max(pick_cand_size,max_degree*_opt_bms_pick_remove);
    }
    pick_cand      = new uint32_t[pick_cand_size];
    
    double runtime = getCpuTime();
    if(_opt_verbosity>0){
        cout<<"c vertex= "<<v_nums<<" ,edges= "<<e_nums<<" ,max_degree= "<<max_degree<<endl;
        cout<<"c --TIME-- build use: "<<runtime<<" s"<<endl;
    }
    last_time_point = runtime;
    // actually we should not use it because readfile and alloc memory use time.
    // and for our paper, FastCDS take the build time in and other competitors donot count this time.
    // this is disadvantage for FastCDS, but FastCDS is good enough even count the time.
    // # but for our open source codes, We don't count that part of the time,
    // because NuCDS, ACO and many other codes do like this.
    // start_time = chrono::steady_clock::now();
    return true;
}

uint32_t FastCDS::fix_points_fast(){
    uint32_t fixed_1_num = 0;
    uint32_t fixed_2_num_clique = 0;
    uint32_t k_clique_degree = max_degree + 1;
    if(_opt_fix_type != 0){
        //method 1 : fix cut point
        uint32_t root               = 1;
        uint32_t root_subtree       = 0;
        uint32_t visit_idx          = 1;
        fill_n(idx_visit,   v_nums+1,   0);
        fill_n(dnf,         v_nums+1,   0);
        my_stack_size   = 1;
        my_stack[0]     = root;
        dnf[root]       = low[root] = visit_idx++;
        uint32_t        fixed_vec_size = 0;
        uint32_t        *fixed_vec;
        fixed_vec       =   new uint32_t[v_nums+1];
        
        while(my_stack_size>0){
            uint32_t top = my_stack[my_stack_size-1];
            if(idx_visit[top] >= degree[top]){
                // cannot to search deeper
                if(--my_stack_size == 0) break;
                uint32_t father = my_stack[my_stack_size-1];
                low[father] = (low[father]<low[top]?low[father]:low[top]);
                if( (father==root && root_subtree > 1) ||
                    (father!=root && low[top]>=dnf[father]) ){
                    fixed_vec[fixed_vec_size++] = father;
                }
            }
            else {
                uint32_t next = adj[top][idx_visit[top]++];
                if(dnf[next]==0){
                    //first visit
                    if(top == root) root_subtree++;
                    dnf[next] = low[next] = visit_idx++;
                    my_stack[my_stack_size++] = next;
                }else{
                    low[top] = (low[top]<dnf[next]?low[top]:dnf[next]);
                }
            }

        }
        for(uint32_t idx=0;idx<fixed_vec_size;++idx){
            uint32_t var = fixed_vec[idx];
            if(fix_type[var] == 0){
                fix_type[var] = 1;
                fix_1_vec[fix_1_vec_size++] = var;
                ++fixed_1_num;
            }
        }
        fill_n(dnf,         v_nums+1,   0);
        delete [] fixed_vec;
    }
    // cout<<"here ok "<<getCpuTime()<<endl;
    // uint32_t ct3 = 0,ct1=0,ct2=0;
    if(_opt_fix_type != 0){
        //method 1 : fix clique points
        uint32_t tmp_map_ct = 0;
        uint32_t *map_ct    = new uint32_t [v_nums+1];
        fill_n(map_ct,v_nums+1,0);
        for(uint32_t var=1;var<=v_nums;++var){
            uint32_t    var_degree      = degree[var];
#ifdef USE_TIME_JUDGE
            if(var%1000==0){
                // cout<<var<<" "<<k_clique_degree<<" "<<ct1<<" "<<ct2<<" "<<ct3<<endl;
                // ct1=ct2=ct3=0;
                if(getCpuTime() > 60){
                    k_clique_degree = 3;
                }
            }
#endif
            if(var_degree < k_clique_degree && fix_type[var]==0 && seen[var]==0){
                if(var_degree == 1){
                    // 2-clique
                    // ct1++;
                    fix_type[var] = 2;
                    fixed_2_num_clique++;
                }else if(var_degree == 2){
                    // 3-clique
                    // ct2++;
                    uint32_t adj1 = adj[var][0];
                    uint32_t adj2 = adj[var][1];
                    if(degree[adj1] < degree[adj2]){
                        for(uint32_t n1_of_adj1 = 0;n1_of_adj1<degree[adj1];n1_of_adj1++){
                            if(adj[adj1][n1_of_adj1] == adj2){
                                fix_type[var] = 2;
                                fixed_2_num_clique ++;

                                seen[adj1] = seen[adj2] = 1;
                                if(degree[adj1] == 2){
                                    fix_type[adj1] = 2;
                                    fixed_2_num_clique ++;
                                }
                                if(degree[adj2] == 2){
                                    fix_type[adj2] = 2;
                                    fixed_2_num_clique ++;
                                }
                                break;
                            }
                        }
                    }else{
                        for(uint32_t n1_of_adj2 = 0;n1_of_adj2<degree[adj2];n1_of_adj2++){
                            if(adj[adj2][n1_of_adj2] == adj1){
                                fix_type[var] = 2;
                                fixed_2_num_clique ++;

                                seen[adj1] = seen[adj2] = 1;
                                if(degree[adj1] == 2){
                                    fix_type[adj1] = 2;
                                    fixed_2_num_clique ++;
                                }
                                if(degree[adj2] == 2){
                                    fix_type[adj2] = 2;
                                    fixed_2_num_clique ++;
                                }
                                break;
                            }
                        }
                    }
                    
                }else{
                    // ct3++;
                    tmp_map_ct++;
                    bool should_fixed_2 = true;
                    for(uint32_t var_idx=0;var_idx<var_degree;++var_idx){
                        uint32_t tmp_var = adj[var][var_idx];
                        if(degree[tmp_var]<var_degree){
                            should_fixed_2 = false;
                            break;
                        }
                        map_ct[tmp_var]=tmp_map_ct;
                    }
                    if(should_fixed_2){
                        for(uint32_t var_idx=0;var_idx<var_degree;++var_idx){
                            uint32_t find_sub_num   = 0;
                            uint32_t should_find_sub= var_degree-var_idx-1;
                            uint32_t tmp_var        = adj[var][var_idx];
                            map_ct[tmp_var]         = 0;
                            uint32_t tmp_var_degree = degree[tmp_var];
                            for(uint32_t n1_idx=0;n1_idx<tmp_var_degree;++n1_idx){
                                uint32_t n1 = adj[tmp_var][n1_idx];
                                if(map_ct[n1] == tmp_map_ct){
                                    if(++find_sub_num >= should_find_sub){
                                        break;
                                    }
                                }
                            }
                            if(find_sub_num < should_find_sub){
                                should_fixed_2 = false;
                                break;
                            }
                        }
                    }
                    if(should_fixed_2){
                        fix_type[var]        = 2;
                        fixed_2_num_clique   ++;
                        uint32_t n1_sz       = degree[var];
                        for(uint32_t n1_idx = 0; n1_idx<n1_sz; ++n1_idx){
                            uint32_t n1     = adj[var][n1_idx];
                            seen[n1]        = 1;
                            if(degree[n1] == var_degree && fix_type[n1]==0){
                                fix_type[n1]    = 2;
                                fixed_2_num_clique++;
                            }
                        }
                    }

                }


            }
        }
        delete [] map_ct;
    }


    fixed_num = fixed_1_num + fixed_2_num_clique;
    double runtime = getCpuTime();
    if(_opt_verbosity > 0){
        cout<<"c Fix type : "<<(uint16_t)_opt_fix_type;
        cout<<" (max_k-clique= "<<k_clique_degree<<")"<<endl;
        cout<<"c Fixed: "<<fixed_num
            <<"("<<fixed_1_num<<","<<fixed_2_num_clique
            <<") out of "<<v_nums<<" vertex"<<endl;
        cout<<"c --TIME-- fix point use: "<<runtime-last_time_point<<" s"<<endl;
    }
    last_time_point = runtime;
    return fixed_num;
}


uint32_t FastCDS::fix_points(){
    uint32_t fixed_1_num = 0;
    uint32_t fixed_2_num_clique = 0;
    uint32_t k_clique_degree = 0;
    if((_opt_fix_type & 0x0001) == (uint8_t)0x0001 || _opt_fix_type == 0x0008){
        //method 1 : fix cut point
        uint32_t root               = 1;
        uint32_t root_subtree       = 0;
        uint32_t visit_idx          = 1;
        fill_n(idx_visit,   v_nums+1,   0);
        fill_n(dnf,         v_nums+1,   0);
        my_stack_size   = 1;
        my_stack[0]     = root;
        dnf[root]       = low[root] = visit_idx++;
        uint32_t        fixed_vec_size = 0;
        uint32_t        *fixed_vec;
        fixed_vec       =   new uint32_t[v_nums+1];
        
        while(my_stack_size>0){
            uint32_t top = my_stack[my_stack_size-1];
            if(idx_visit[top] >= degree[top]){
                // cannot to search deeper
                if(--my_stack_size == 0) break;
                uint32_t father = my_stack[my_stack_size-1];
                low[father] = (low[father]<low[top]?low[father]:low[top]);
                if( (father==root && root_subtree > 1) ||
                    (father!=root && low[top]>=dnf[father]) ){
                    fixed_vec[fixed_vec_size++] = father;
                }
            }
            else {
                uint32_t next = adj[top][idx_visit[top]++];
                if(dnf[next]==0){
                    //first visit
                    if(top == root) root_subtree++;
                    dnf[next] = low[next] = visit_idx++;
                    my_stack[my_stack_size++] = next;
                }else{
                    low[top] = (low[top]<dnf[next]?low[top]:dnf[next]);
                }
            }

        }
        for(uint32_t idx=0;idx<fixed_vec_size;++idx){
            uint32_t var = fixed_vec[idx];
            if(fix_type[var] == 0){
                fix_type[var] = 1;
                fix_1_vec[fix_1_vec_size++] = var;
                ++fixed_1_num;
            }
        }
        fill_n(dnf,         v_nums+1,   0);
        delete [] fixed_vec;
    }
    // cout<<"here ok "<<getCpuTime()<<endl;
    // uint32_t ct3 = 0,ct1=0,ct2=0;
    // nodes with degree=k in a k-clique should be fixed-0,
    // unless all other have fixed-0, but it is not possible to be a complete graph for CDS instances.
    if( ( (_opt_fix_type & 0x000e) != 0x0000 ) &&\
        (max_degree<1.5*1000*1000 || avg_degree<18) ){
        
        // calculate the proper k in k-clique.
        if((_opt_fix_type & 0x0008) == 0x0008){
            k_clique_degree = (uint16_t)_opt_fix_type>>4;
            if(k_clique_degree == 0){
                // if construct time is bigger than 300s, shutdown.
                // use conditions replace if(runtime>300s)
                if(    (e_nums      > 100*1000*1000)
                    || (max_degree  > 200*1000)
                    || (v_nums      > 58*1000*1000)
                    || (v_nums>600*1000  && max_degree>84*1000 && avg_degree>20.0)
                    || (v_nums>420*1000  && max_degree>45*1000 && avg_degree>40.0)
                    || (v_nums>850*1000  && max_degree>65*1000 && avg_degree>37.0)
                    || (v_nums>1300*1000 && max_degree>20*1000 && avg_degree>19.0)
                    || (v_nums>2100*1000 && max_degree>97*1000 && avg_degree>16.0)
                    || (v_nums>1600*1000 && max_degree>35*1000 && avg_degree>13.0)
                ){
                    k_clique_degree = 3;
                }else{
                    k_clique_degree = max_degree+1;
                }
            }
        }else if((_opt_fix_type & 0x0004) == 0x0004){
            k_clique_degree = max_degree + 1;
        }else if((_opt_fix_type & 0x0002) == 0x0002){
            k_clique_degree = 3;
        }
#ifdef USE_TIME_JUDGE
        k_clique_degree = max_degree + 1;
#endif
        // uint8_t *seen = new uint8_t[v_nums+1];
        fill_n(seen,v_nums+1,0);
        for(uint32_t var=1;var<=v_nums;++var){
            uint32_t    var_degree      = degree[var];
#ifdef USE_TIME_JUDGE
            if(var%1000==0){
                // cout<<var<<" "<<k_clique_degree<<" "<<ct1<<" "<<ct2<<" "<<ct3<<endl;
                // ct1=ct2=ct3=0;
                if(getCpuTime() > 60){
                    k_clique_degree = 3;
                }
            }
#endif
            if(var_degree<k_clique_degree && fix_type[var] == 0 && seen[var] == 0){
                
                if(var_degree == 1){
                    // 2-clique
                    // ct1++;
                    fix_type[var]    = 2;
                    fixed_2_num_clique++;
                    continue;
                }
                else if(var_degree == 2){
                    // 3-clique
                    // ct2++;
                    uint32_t adj1 = adj[var][0];
                    uint32_t adj2 = adj[var][1];
                    for(uint32_t n1_of_adj1 = 0;n1_of_adj1<degree[adj1];n1_of_adj1++){
                        if(adj[adj1][n1_of_adj1] == adj2){
                            fix_type[var] = 2;
                            fixed_2_num_clique ++;
                            seen[adj1] = seen[adj2] = 1;
                            if(degree[adj1] == 2){
                                fix_type[adj1] = 2;
                                fixed_2_num_clique ++;
                            }
                            if(degree[adj2] == 2){
                                fix_type[adj2] = 2;
                                fixed_2_num_clique ++;
                            }
                            break;
                        }
                    }
                    continue;
                }
        
        
                else{
                    // ct3++;
                    bool        should_fixed_2  = true;
                    
                    // k-clique with k>3
                    for(uint32_t n1_idx = 0;n1_idx<var_degree-1;++n1_idx){
                        uint32_t n1 = adj[var][n1_idx];
                        if(degree[n1] < var_degree){
                            should_fixed_2 = false;
                            goto break_point;
                        }else{
                            unordered_set<uint32_t> is_in;
                            for(uint32_t n1_idx_inner = 0; n1_idx_inner<degree[n1]; ++n1_idx_inner){
                                is_in.insert(adj[n1][n1_idx_inner]);
                            }
                            for(uint32_t n1_idx_inner = n1_idx+1; n1_idx_inner<var_degree; ++n1_idx_inner){
                                uint32_t n1_inner = adj[var][n1_idx_inner];
                                if(is_in.find(n1_inner) == is_in.end()){
                                    should_fixed_2 = false;
                                    goto break_point;
                                }
                            }
                        }
                    }
                    break_point:;
                    if(should_fixed_2){
                        //fixed self and all k-degree nodes
                        fix_type[var]        = 2;
                        fixed_2_num_clique   ++;
                        for(uint32_t n1_idx = 0; n1_idx<degree[var]; ++n1_idx){
                            uint32_t n1     = adj[var][n1_idx];
                            seen[n1]        = 1;
                            if(degree[n1] == var_degree && fix_type[n1]==0){
                                fix_type[n1]    = 2;
                                fixed_2_num_clique++;
                            }
                        }
                    }
                }
            }
        }
        // delete [] seen;
        
    }

    
fix_end:;
    fixed_num = fixed_1_num + fixed_2_num_clique;
    double runtime = getCpuTime();
    if(_opt_verbosity > 0){
        cout<<"c Fix type : "<<(uint16_t)_opt_fix_type;
        cout<<" (max_k-clique= "<<k_clique_degree<<")"<<endl;
        cout<<"c Fixed: "<<fixed_num
            <<"("<<fixed_1_num<<","<<fixed_2_num_clique
            <<") out of "<<v_nums<<" vertex"<<endl;
        cout<<"c --TIME-- fix point use: "<<runtime-last_time_point<<" s"<<endl;
    }
    last_time_point = runtime;
    return fixed_num;
}

void FastCDS::renew_cut(){
    // used before every pick_remove_SUB;
    for(uint32_t i=0; i<cut_vec_size; ++i) is_cut[cut_vec[i]]=0;
    cut_vec_size = 0;
    for(uint32_t i=0;i<cds_size;++i) idx_visit[cds_vec[i]]=0;

    fill_n(dnf,v_nums+1,0);
    uint32_t r = rand()%cds_size;
    uint32_t    root            =   cds_vec[r];
    uint32_t    root_subtree    =   0;
    uint32_t    visit_idx       =   1;
    
    my_stack_size   = 1;
    my_stack[0]     = root;
    dnf[root]       = low[root] = visit_idx++;
    
    while(my_stack_size>0){
        uint32_t top = my_stack[my_stack_size-1];
        if(idx_visit[top] >= degree[top]){
            // cannot to search deeper
            if(--my_stack_size == 0) break;
            uint32_t father = my_stack[my_stack_size-1];
            low[father] = (low[father]<low[top]?low[father]:low[top]);
            if( (father==root && root_subtree > 1) ||
                (father!=root && low[top]>=dnf[father]) ){
                is_cut[father] = 1;
                cut_vec[cut_vec_size++] = father;
            }
        }
        else {
            uint32_t next = adj[top][idx_visit[top]++];
            if(dnf[next]==0){
                //first visit
                if(cds[next]==0) continue;
                if(top == root) root_subtree++;
                dnf[next] = low[next] = visit_idx++;
                my_stack[my_stack_size++] = next;
            }else{
                low[top] = (low[top]<dnf[next]?low[top]:dnf[next]);
            }
        }

    }
    
}


void FastCDS::construct_tree_random_based(){
    for(uint32_t idx=0;idx<leaf_vec_size;++idx){is_leaf[leaf_vec[idx]]=0;}
    
    fill_n(seen,v_nums+1,0);
    fill_n(son_num,v_nums+1,0);
    leaf_vec_size   = 0;
    
    q_front          = 0;
    q_back           = 1;
    // rand pick k times, to get the max degree one;
    root             = cds_vec[rand()%cds_size];
    // for(uint32_t ct=1;ct<3;++ct){
    //     uint32_t rand_pick  = cds_vec[rand()%cds_size];
    //     if(degree[rand_pick]>degree[root]){
    //         root = rand_pick;
    //     }
    // }
    
    my_queue[0]            = root;
    seen[root]             = 1;
    while (q_front < q_back) {
        // uint32_t front_node     = my_queue[q_front++];
        //diversification
        uint32_t pick_node_idx  = q_front+rand()%(q_back-q_front);
        uint32_t front_node     = my_queue[pick_node_idx];
        my_queue[pick_node_idx] = my_queue[q_front++];

        // pick son is not randomlized. could be improved.
        for(uint32_t n1_idx=0;n1_idx<degree[front_node];++n1_idx){
            uint32_t n1 = adj[front_node][n1_idx];
            if(cds[n1]==1 && seen[n1]==0){
                son_num[front_node]++;
                father[n1] = front_node;
                my_queue[q_back++] = n1;
                seen[n1]    = 1;
            }
        }
        if(son_num[front_node] == 0){
            is_leaf[front_node]         = 1;
            idx_in_leaf_vec[front_node] = leaf_vec_size;
            leaf_vec[leaf_vec_size++]   = front_node;
        }
    }
    
   
}




void FastCDS::construct_tree_random_based_for_big(){
    for(uint32_t idx=0;idx<leaf_vec_size;++idx){is_leaf[leaf_vec[idx]]=0;}
    
    fill_n(seen,v_nums+1,0);
    fill_n(son_num,v_nums+1,0);
    leaf_vec_size   = 0;
    
    q_front          = 0;
    q_back           = 1;
    root             = cds_vec[rand()%cds_size];
    
    my_queue[0]            = root;
    seen[root]             = 1;
    while (q_front < q_back) {
        // randomly pick next node.
        uint32_t pick_node_idx  = q_front+rand()%(q_back-q_front);
        uint32_t front_node     = my_queue[pick_node_idx];
        my_queue[pick_node_idx] = my_queue[q_front++];
        seen_before_queue[front_node] = 1;


        // random pick a father for this node.
        best_father_cand_size = 0;
        uint32_t best_father_degree = 0;
        uint32_t n1_neighbor_sz = degree[front_node];
        for(uint32_t n1_idx=0;n1_idx<n1_neighbor_sz;++n1_idx){
            uint32_t n1 = adj[front_node][n1_idx];
            if(cds[n1]==1){
                if(seen_before_queue[n1]==1){
                    if(degree[n1] > best_father_degree){
                        best_father_cand_size = 1;
                        best_father_cand[0] = n1;
                    }else if(degree[n1] == best_father_degree)
                        best_father_cand[best_father_cand_size++] = n1;
                }
                if(seen[n1]==0){
                    my_queue[q_back++] = n1;
                    seen[n1] = 1;
                }
            }
        }
        if(front_node == root){continue;}
        uint32_t rand_pick_father = best_father_cand[rand()%best_father_cand_size];
        son_num[rand_pick_father]++;
        father[front_node] = rand_pick_father;
    }


    for(uint32_t cds_idx=0;cds_idx<cds_size;++cds_idx){
        uint32_t cds_node = cds_vec[cds_idx];
        seen_before_queue[cds_node] = 0;
        if(son_num[cds_node] == 0){
            is_leaf[cds_node] = 1;
            idx_in_leaf_vec[cds_node] = leaf_vec_size;
            leaf_vec[leaf_vec_size++] = cds_node;
        }
    }
   
}



void  FastCDS::set_subgraph_data_structure(){
    uint32_t v_mems = v_nums+1;
    subgraph_score  = new uint32_t[v_mems];
    subgraph_adj    = new uint32_t*[v_mems];
    subgraph_degree = new uint32_t[v_mems];
    subgraph_is_grey= new uint8_t[v_mems];
    subgraph_dom_num= new uint32_t[v_mems];
    subgraph_dom_by = new uint32_t[v_mems];
    for(uint32_t v_idx=1;v_idx<=v_nums;++v_idx){
        subgraph_adj[v_idx]= new uint32_t[degree[v_idx]+1];
    }
    subgraph_heap = new My_heap<int32_t>(v_nums,[&](uint32_t a, uint32_t b)->bool{
        if(subgraph_score[a] < subgraph_score[b]){
            return true;
        }else if(subgraph_score[a] == subgraph_score[b]){
            if(freq_score[a] > freq_score[b]){
                return true;
            }else if (freq_score[a]==freq_score[b] && pick_swap_num[a]<pick_swap_num[b]){
                return true;
            }
        }
        return false;
    });
    fill_n(subgraph_dom_num,v_mems,0);
    fill_n(subgraph_is_grey,v_mems,0);
    fill_n(subgraph_degree,v_mems,0);
}

void FastCDS::releases_subgraph_data_structure(){
    delete subgraph_heap;
    delete [] subgraph_degree;
    delete [] subgraph_score;
    delete [] subgraph_is_grey;
    delete [] subgraph_dom_num;
    delete [] subgraph_dom_by;
    for(uint32_t v_idx;v_idx<=v_nums;v_idx++){
        delete [] subgraph_adj[v_idx];
    }
    delete [] subgraph_adj;
}


void FastCDS::construct_tree_subgraph_based(){
    for(uint32_t idx=0;idx<leaf_vec_size;++idx){is_leaf[leaf_vec[idx]]=0;}
    fill_n(seen,v_nums+1,0); // a node who has choosen in tree
    fill_n(son_num,v_nums+1,0);
    
    leaf_vec_size   = 0;
    subgraph_heap->my_heap_size = 0;
    uint32_t max_degree_subvertex = 0;
    int32_t  max_subgraph_degree = 0;
    for(uint32_t sub_idx_v=0;sub_idx_v<cds_size;++sub_idx_v){
        uint32_t sub_vertex = cds_vec[sub_idx_v];
        subgraph_degree[sub_vertex] = 0;
        uint32_t adj_sz = degree[sub_vertex];
        uint32_t adj_cds_sz = 0;
        for(uint32_t adj_idx=0;adj_idx<adj_sz;++adj_idx){
            uint32_t adj_vertex = adj[sub_vertex][adj_idx];
            if(cds[adj_vertex] == 1){
                subgraph_adj[sub_vertex][adj_cds_sz++] = adj_vertex;
            }
        }
        subgraph_degree[sub_vertex] = adj_cds_sz;
        subgraph_score[sub_vertex]  = adj_cds_sz + 1;
        if(adj_cds_sz > max_subgraph_degree){
            max_subgraph_degree = adj_cds_sz;
            max_degree_subvertex = sub_vertex;
        }
    }

    // pick the max_degree node as root;
    // root   = max_degree_subvertex;
    // uint32_t root_degree = max_subgraph_degree;

    // random pick root
    root    = cds_vec[rand()%cds_size];
    uint32_t root_degree = subgraph_degree[root];

    subgraph_score[root]--;
    subgraph_is_grey[root] = 1;
    for(uint32_t n1_idx=0;n1_idx<root_degree;++n1_idx){
        subgraph_score[subgraph_adj[root][n1_idx]]--;
    }
    
    subgraph_heap->add_node(root);
    uint32_t    add_node    =   0;
    uint32_t    sub_dom_sz  =   1;
    while(sub_dom_sz < cds_size){
        add_node = subgraph_heap->pop_top();
        while(subgraph_score[add_node]==0) add_node = subgraph_heap->pop_top();
        seen[add_node] = 1;
        if(add_node != root){
            // deal with leaves (most of the subCDS nodes should not be leaves)
            uint32_t max_degree_of_father = 0;
            best_father_cand_size = 0;
            uint32_t father_search_degree = subgraph_degree[add_node];
            for(uint32_t idx_father=0;idx_father<father_search_degree;++idx_father){
                uint32_t search_father = subgraph_adj[add_node][idx_father];
                if(seen[search_father]){
                    uint32_t search_degree = subgraph_degree[search_father]; 
                    if(search_degree>max_degree_of_father){
                        max_degree_of_father = search_degree;
                        best_father_cand_size = 1;
                        best_father_cand[0] = search_father;
                    }else if(search_degree == max_degree_of_father){
                        best_father_cand[best_father_cand_size++] = search_father;
                    }
                }
            }
            uint32_t rand_picked_father = best_father_cand[rand()%best_father_cand_size];
            son_num[rand_picked_father]++;
            father[add_node] = rand_picked_father;
        }
        sub_dom_sz += subgraph_score[add_node];
        subgraph_is_grey[add_node] = 0;
        uint32_t v_score = subgraph_score[add_node];
        if(++subgraph_dom_num[add_node]==2) subgraph_score[subgraph_dom_by[add_node]]++;
        else if(subgraph_dom_num[add_node]==1) subgraph_dom_by[add_node]=add_node;
        uint32_t add_node_degree = subgraph_degree[add_node];
        for(uint32_t n1_idx=0;n1_idx<add_node_degree;++n1_idx){
            uint32_t n1 = subgraph_adj[add_node][n1_idx];
            subgraph_dom_num[n1]++;
            if(subgraph_dom_num[n1]==1){
                subgraph_is_grey[n1] = 1;
                subgraph_dom_by[n1]  = add_node;
                subgraph_score[n1]--;
                subgraph_heap->add_node(n1);
                uint32_t n1_degree = subgraph_degree[n1];
                for(uint32_t n2_idx=0;n2_idx<n1_degree;++n2_idx){
                    uint32_t n2 = subgraph_adj[n1][n2_idx];
                    if(seen[n2]==0){
                        subgraph_score[n2]--;
                        if(subgraph_is_grey[n2]){
                            subgraph_heap->conduct_down(subgraph_heap->idx_in_my_heap[n2]);
                        }
                    }
                }
            }else if(subgraph_dom_num[n1]==2){
                subgraph_score[subgraph_dom_by[n1]]++;
            }
        }
        if(add_node != root){
            subgraph_score[add_node] = -v_score;
        }else{
            subgraph_score[add_node] = -v_score-1;
        }
    }

    for(uint32_t idx_subgraph=0;idx_subgraph<cds_size;++idx_subgraph){
        uint32_t subgraph_node = cds_vec[idx_subgraph];
        if(seen[subgraph_node]==0){
            //pick father for rest leaf nodes 
            uint32_t max_degree_of_father = 0;
            best_father_cand_size = 0;
            uint32_t father_search_degree = subgraph_degree[subgraph_node];
            for(uint32_t idx_father=0;idx_father<father_search_degree;++idx_father){
                uint32_t search_father = subgraph_adj[subgraph_node][idx_father];
                if(seen[search_father]){
                    uint32_t search_degree = subgraph_degree[search_father]; 
                    if(search_degree>max_degree_of_father){
                        max_degree_of_father = search_degree;
                        best_father_cand_size = 1;
                        best_father_cand[0] = search_father;
                    }else if(search_degree == max_degree_of_father){
                        best_father_cand[best_father_cand_size++] = search_father;
                    }
                }
            }
            uint32_t rand_picked_father = best_father_cand[rand()%best_father_cand_size];
            son_num[rand_picked_father]++;
            father[subgraph_node] = rand_picked_father;
        }
    }
    for(uint32_t idx_subgraph=0;idx_subgraph<cds_size;++idx_subgraph){
        uint32_t subgraph_node = cds_vec[idx_subgraph];
        if(son_num[subgraph_node] == 0){
            is_leaf[subgraph_node] = 1;
            idx_in_leaf_vec[subgraph_node] = leaf_vec_size;
            leaf_vec[leaf_vec_size++] = subgraph_node;
        }
    }


    for(uint32_t sub_idx_v=0;sub_idx_v<cds_size;++sub_idx_v){
        uint32_t sub_vertex = cds_vec[sub_idx_v];
        subgraph_is_grey[sub_vertex] = 0;
        subgraph_degree[sub_vertex] = 0;
        subgraph_dom_num[sub_vertex] = 0;
    }

}


void FastCDS::construct_tree_hueristic_based(){
    for(uint32_t idx=0;idx<leaf_vec_size;++idx){is_leaf[leaf_vec[idx]]=0;}
    
    fill_n(seen,v_nums+1,0);
    fill_n(son_num,v_nums+1,0);
    leaf_vec_size   = 0;
    
    q_front          = 0;
    q_back           = 1;
    // rand pick k times, to get the max degree one;
    root             = cds_vec[rand()%cds_size];
    // for(uint32_t ct=1;ct<3;++ct){
    //     uint32_t rand_pick  = cds_vec[rand()%cds_size];
    //     if(degree[rand_pick]>degree[root]){
    //         root = rand_pick;
    //     }
    // }
    seen[root]             = 1;

    My_heap<int32_t> my_heap(v_nums,[&](uint32_t a, uint32_t b)->bool{
        if(score[a] < score[b]){
            return true;
        }else if(score[a] == score[b]){
            if(freq_score[a] > freq_score[b]){
                return true;
            }else if (freq_score[a]==freq_score[b] && pick_swap_num[a]<pick_swap_num[b]){
                return true;
            }
        }
        return false;
    });
    my_heap.add_node(root);
    while(my_heap.my_heap_size>0){
        uint32_t pick_node      = my_heap.pop_top();
        uint32_t pick_degree    = degree[pick_node];
        for(uint32_t n1_idx=0;n1_idx<pick_degree;++n1_idx){
            uint32_t n1 = adj[pick_node][n1_idx];
            if(cds[n1]==1 && seen[n1]==0){
                son_num[pick_node]++;
                father[n1] = pick_node;
                my_heap.add_node(n1);
                seen[n1]   = 1;
            }
        }
        if(son_num[pick_node] == 0){
            is_leaf[pick_node]         = 1;
            idx_in_leaf_vec[pick_node] = leaf_vec_size;
            leaf_vec[leaf_vec_size++]  = pick_node;
        }
    }
    
    
   
}



// The initial soln using bms is really bad.
// Do not use it !!!!
uint32_t FastCDS::construct_1hop_bms(){
    cout<<"c This is a bad idea"<<endl;
    exit(20);
    /*
    fix_points();
        
    fill_n(cds,         v_nums+1,   0);
    fill_n(is_grey,     v_nums+1,   0);
    is_grey[max_degree_vertex] = 1;
    
    for(uint32_t idx=1;idx<=v_nums;++idx) score[idx] = degree[idx]+1;
    score[max_degree_vertex]--;
    for(uint32_t n1_idx=0;n1_idx<degree[max_degree_vertex];++n1_idx) score[adj[max_degree_vertex][n1_idx]]--;
    
    // only save grey nodes whose score > 0
    // no need to deal with idx because idx is rand selected and can be saved during BMS.
    my_stack_size = 0;
    my_stack[my_stack_size++] = max_degree_vertex;
    

    
    int32_t     add_node            = -1;
    uint32_t    add_node_stack_idx  = 0;
    uint32_t    dominated_sz        = 1;
                cds_size            = 0;
                white_vec_size      = 0;
                grey_vec_size       = 0;
    uint32_t    wait_add_list_size  = 0;
    uint32_t    *wait_add_list      = new uint32_t[fixed_num];
    while(dominated_sz < v_nums){
        
        bool from_add_list = false;
        if(wait_add_list_size>0){
            add_node = wait_add_list[--wait_add_list_size];
            from_add_list = true;
        }else{
            // use BMS to select the add_node;
            // the bigger bms_select_num, the better.
            int32_t max_score = -1;
            if(my_stack_size  <= _opt_bms_construct * 1.0){
                // visit all nodes;
                for(uint16_t idx=0; idx<my_stack_size; ++idx){
                    uint32_t tmp_node = my_stack[idx];
                    if(score[tmp_node] > max_score){
                        max_score           = score[tmp_node];
                        add_node            = tmp_node;
                        add_node_stack_idx  = idx;
                    }
                }
            }else{
                //use BMS
                for(uint16_t bms_idx=0; bms_idx<_opt_bms_construct; ++bms_idx){
                    uint32_t  idx      = rand()%my_stack_size;
                    uint32_t tmp_node  = my_stack[idx];
                    if(score[tmp_node] > max_score){
                        max_score           = score[tmp_node];
                        add_node            = tmp_node;
                        add_node_stack_idx  = idx;
                    }
                }
            }
        }
        
        if(add_node > 0){
            cds_vec[cds_size]   = add_node;
            idx_in_cds_vec[add_node]  = cds_size++;
            dominated_sz        += score[add_node];
            is_grey[add_node]   = 0;
            cds[add_node]       = 1;
            if(!from_add_list){
                my_stack[add_node_stack_idx] = my_stack[--my_stack_size];
            }
            uint32_t v_score    = score[add_node];
            if(++dom_num[add_node] == 2) score[dom_by[add_node]]++;
            else if(dom_num[add_node]==1) dom_by[add_node] = add_node;
            // update n2 score of n1 (n1 white -> grey)
            for(uint32_t n1_idx=0;n1_idx<degree[add_node];++n1_idx){
                uint32_t n1 = adj[add_node][n1_idx];
                dom_num[n1]++;
                if(dom_num[n1]==1){
                    //white -> grey
                    is_grey[n1] = 1;
                    dom_by[n1]  = add_node;
                    if(fix_type[n1] == 1){
                        wait_add_list[wait_add_list_size++] = n1;
                    }else if(fix_type[n1] == 0){
                        my_stack[my_stack_size++] = n1;
                    }
                    score[n1]--;
                    for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                        uint32_t n2 = adj[n1][n2_idx];
                        if(cds[n2] == 0){
                            score[n2]--;
                        }
                    }
                }else if(dom_num[n1]==2){
                    score[dom_by[n1]]++;
                }
            }
            if(add_node != max_degree_vertex){
                score[add_node] = -v_score;
            }else{
                score[add_node] = -v_score-1;
            }
        }
    }
    for(uint32_t vertex=1;vertex<=v_nums;++vertex){
        if(cds[vertex]==0){
            grey_vec[grey_vec_size] = vertex;
            idx_in_grey_vec[vertex] = grey_vec_size++;
        }
    }
    undom_size = 0;
    delete [] wait_add_list;
    copy(score, score+v_nums+1, freq_score);
    update_best_soln();
    if(_opt_verbosity>0){
        cout<<"c Initial cds_size = "<<best_cds_size<<endl;
        double runtime = getCpuTime();
        cout<<"c --TIME-- construct(BMS,"<<_opt_bms_construct<<") use: "<<runtime-last_time_point<<" s"<<endl;
        last_time_point = runtime;
    }
    return cds_size;
    */
}


uint32_t FastCDS::construct_1hop(){
    // fix_points();
    fix_points_fast();
    fill_n(cds,         v_nums+1,   0);
    fill_n(is_grey,     v_nums+1,   0);

    is_grey[max_degree_vertex] = 1;
    for(uint32_t idx=1;idx<=v_nums;++idx) score[idx] = degree[idx]+1;
    score[max_degree_vertex]--;
    for(uint32_t n1_idx=0;n1_idx<degree[max_degree_vertex];++n1_idx) score[adj[max_degree_vertex][n1_idx]]--;
#ifdef LOOKAHEAD_PICK_1HOP
    uint32_t        look_ahead_list_size = 0;
    const uint32_t  look_ahead_tenuer    = 10;
    uint32_t        *look_ahead_list     =  new uint32_t[look_ahead_tenuer+10];
#endif

    // only save grey nodes whose score > 0 and fix_type!=2
#ifdef USE_RAND_IN_CONSTRUCTION
    My_heap<uint32_t> heap_degree(v_nums,degree);
    heap_degree.add_node(max_degree_vertex);
    int8_t      mod                 = 1;
    uint32_t    *new_cover_nodes    = new uint32_t[v_nums+1];
    uint32_t    new_cover_nodes_sz  = 0;
#endif
    My_heap<int32_t>  heap_score (v_nums,score);
    heap_score.add_node(max_degree_vertex);
    
    int32_t     add_node            = -1;
    uint32_t    dominated_sz        = 1;
                cds_size            = 0;
                white_vec_size      = 0;
                grey_vec_size       = 0;
    uint32_t    wait_add_list_size  = 0;
    uint32_t    *wait_add_list      = new uint32_t[fixed_num];
    while(dominated_sz < v_nums){
        
        bool from_add_list = false;
        if(wait_add_list_size>0){
            add_node = wait_add_list[--wait_add_list_size];
            // from_add_list = true;
        }else{
#ifdef USE_RAND_IN_CONSTRUCTION
            if(mod == 1){
#endif

#ifdef LOOKAHEAD_PICK_1HOP
                // pick the top node randomly use look-ahead
                look_ahead_list_size = 0;
                uint32_t visited_idx = 0;
                uint32_t top_var     = heap_score.my_heap[0];
                uint32_t top_score   = score[top_var];
                look_ahead_list[look_ahead_list_size++] = top_var;
                while(visited_idx < look_ahead_list_size){
                    uint32_t pick_node_list = look_ahead_list[visited_idx];
                    uint32_t idx_pick_node  = heap_score.idx_in_my_heap[pick_node_list];
                    uint32_t son_idx = heap_score.left_son(idx_pick_node);
                    if(son_idx >= heap_score.my_heap_size){break;}
                    uint32_t son     = heap_score.my_heap[son_idx];
                    if(score[son] == top_score) look_ahead_list[look_ahead_list_size++] = son;
                    son_idx = heap_score.right_son(idx_pick_node);
                    if(son_idx >= heap_score.my_heap_size){break;}
                    son     = heap_score.my_heap[son_idx];
                    if(score[son] == top_score) look_ahead_list[look_ahead_list_size++] = son;
                    if(look_ahead_list_size>look_ahead_tenuer) break;
                    visited_idx++;
                }
                add_node    = look_ahead_list[rand()%look_ahead_list_size];
                heap_score.remove_node(add_node);
#ifdef USE_RAND_IN_CONSTRUCTION
                heap_degree.remove_node(add_node);
#endif               
                
#else
                // pick the top node
                add_node    = heap_score.my_heap[0];
                heap_score. remove_node(add_node);
#ifdef USE_RAND_IN_CONSTRUCTION
                heap_degree.remove_node(add_node);
#endif
#endif

#ifdef USE_RAND_IN_CONSTRUCTION
                new_cover_nodes_sz = 0;
            }else if(true || new_cover_nodes_sz>0){
                // for some possibility, add a grey node.
                double harmony = 0;
                for(uint32_t i=1;i<=degree[heap_degree.my_heap[0]];++i) harmony += 1/(i+0.0);
                harmony = pow(harmony,-0.5);
                // find the rand maybe not so useful for big instance. the less possibility, the better result.
                if((double)(rand()%INT16_MAX)/INT16_MAX < harmony){
                    // rand pick a grey node in newly coverd nodes
                    add_node = new_cover_nodes[rand()%new_cover_nodes_sz];
                    // rand pick a grey node with score > 0;
                    // add_node = heap_score.my_heap[rand() % heap_score.my_heap_size];
                    // while(score[add_node] < 1){
                    //     heap_score. remove_node(add_node);
                    //     heap_degree.remove_node(add_node);
                    //     add_node = heap_score.my_heap[rand() % heap_score.my_heap_size];
                    // }
                }else add_node = -1;
            }
            mod = -mod;
#endif
        }

        
        if(add_node > 0){
            cds_vec[cds_size]   = add_node;
            idx_in_cds_vec[add_node]  = cds_size++;
            dominated_sz        += score[add_node];
            is_grey[add_node]   = 0;
            cds[add_node]       = 1;
            uint32_t v_score    = score[add_node];
            if(++dom_num[add_node]==2) score[dom_by[add_node]]++;
            else if(dom_num[add_node]==1) dom_by[add_node] = add_node;
            // update n2 score of n1 (n1 white -> grey)
            uint32_t add_node_degree = degree[add_node];
            for(uint32_t n1_idx=0; n1_idx<add_node_degree; ++n1_idx){
                uint32_t n1 = adj[add_node][n1_idx];
                dom_num[n1]++;
                if(dom_num[n1]==1){
                    //white -> grey
#ifdef USE_RAND_IN_CONSTRUCTION
                    new_cover_nodes[new_cover_nodes_sz++] = n1;    
#endif
                    is_grey[n1] = 1;
                    dom_by[n1]  = add_node;
                    score[n1]--;
                    if(fix_type[n1] == 1){
                        wait_add_list[wait_add_list_size++] = n1;
                    }else if(fix_type[n1] == 0){
                        heap_score. add_node(n1);
#ifdef USE_RAND_IN_CONSTRUCTION
                        heap_degree.add_node(n1);
#endif
                    }
                    for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                        uint32_t n2 = adj[n1][n2_idx];
                        if(cds[n2]==0){
                            score[n2]--;
                            if(fix_type[n2]!=2 && is_grey[n2])
                                heap_score.conduct_down(heap_score.idx_in_my_heap[n2]);
                        }
                    }
                }else if(dom_num[n1]==2){
                    score[dom_by[n1]]++;
                }
            }
            if(add_node != max_degree_vertex){
                score[add_node] = -v_score;
            }else{
                score[add_node] = -v_score-1;
            }
        }
    }
    for(uint32_t vertex=1;vertex<=v_nums;++vertex){
        if(cds[vertex]==0){
            grey_vec[grey_vec_size] = vertex;
            idx_in_grey_vec[vertex] = grey_vec_size++;
        }
    }
    undom_size = 0;
    delete [] wait_add_list;
#ifdef USE_RAND_IN_CONSTRUCTION
    delete [] new_cover_nodes;
#endif
#ifdef LOOKAHEAD_PICK_1HOP
    delete [] look_ahead_list;
#endif
    copy(score, score+v_nums+1, freq_score);
    update_best_soln();
    double runtime = getCpuTime();
    if(_opt_verbosity>0){
        cout<<"c Initial cds_size = "<<best_cds_size<<endl;
        cout<<"c --TIME-- construct(HEAP) use: "<<runtime-last_time_point<<" s"<<endl;
    }
    last_time_point = runtime;
    if(_opt_only_output_init_soln!=0){
        cout<<"v ";
        for(uint32_t i=1;i<=v_nums;++i){
            if(cds[i] == 1) cout<<i<<" ";
        }cout<<endl;
        cout<<"c Print Init Soln and Exit!"<<endl;
        exit(20);
    }
    return cds_size;
    
}


void FastCDS::reconstruct_with_upper_bound(uint32_t upper_bound){
    fill_n(cds,         v_nums+1,   0);
    fill_n(is_grey,     v_nums+1,   0);
    // fill_n(tabu_remove, v_nums+1,   0);
    fill_n(conf_change, v_nums+1,   1);
    fill_n(dom_num,     v_nums+1,   0);
    //pick root and initialize the score.
    uint32_t rand_pick_root = rand()%v_nums+1;
    // for(uint32_t i=0;i<10;++i){
    //     uint32_t tmp_pick = rand()%v_nums;
    //     if(degree[rand_pick_root]<degree[tmp_pick] && fix_type[tmp_pick]!=2){
    //         rand_pick_root = tmp_pick;
    //     }
    // }
    is_grey[rand_pick_root] = 1;
    for(uint32_t idx=1;idx<=v_nums;++idx) {
        score[idx] = degree[idx]+1;
        freq_score[idx] = freq[idx];
        for(uint32_t n1_idx=0;n1_idx<degree[idx];++n1_idx){
            freq_score[idx] += freq[adj[idx][n1_idx]];
        }
    }
    score[rand_pick_root]--;
    freq_score[rand_pick_root]-=freq[rand_pick_root];
    for(uint32_t n1_idx=0;n1_idx<degree[rand_pick_root];++n1_idx){
        uint32_t n1 = adj[rand_pick_root][n1_idx];
        score[n1]--;
        freq_score[n1]-=freq[rand_pick_root];
    }
    //only save grey nodes whose freq_score>0 and fix_type!=2
    My_heap<int32_t> heap_fscore(v_nums,freq_score);
    heap_fscore.add_node(rand_pick_root);
    
    int32_t     add_node        = -1;
    uint32_t    dominated_sz    = 1;
                white_vec_size  = 0;
                cds_size        = 0;
                grey_vec_size   = 0;
    uint32_t    wait_add_list_size  = 0;
    uint32_t    *wait_add_list      = new uint32_t[fixed_num];
    
    while(dominated_sz < v_nums && cds_size<upper_bound){
            
            bool from_add_list = false;
            if(wait_add_list_size>0){
                add_node = wait_add_list[--wait_add_list_size];
                from_add_list = true;
            }else{
                add_node      = heap_fscore.my_heap[0];
            }
            
            
            cds_vec[cds_size]   = add_node;
            idx_in_cds_vec[add_node]  = cds_size++;
            dominated_sz        += score[add_node];
            is_grey[add_node]   = 0;
            cds[add_node]       = 1;
            if(!from_add_list){
                heap_fscore. remove_node(add_node);
            }
            uint32_t v_score    = score[add_node];
            uint32_t v_fscore   = freq_score[add_node];
            if(++dom_num[add_node]==2) {
                uint32_t v_dom_by = dom_by[add_node];
                ++score[v_dom_by];
                freq_score[v_dom_by] += freq[add_node];
            }
            else if(dom_num[add_node]==1) dom_by[add_node] = add_node;
            // update n2 score of n1 (n1 white -> grey)
            for(uint32_t n1_idx=0;n1_idx<degree[add_node];++n1_idx){
                uint32_t n1 = adj[add_node][n1_idx];
                dom_num[n1]++;
                if(dom_num[n1]==1){
                    //white -> grey
                    is_grey[n1] = 1;
                    dom_by[n1]  = add_node;
                    if(fix_type[n1] == 1){
                        wait_add_list[wait_add_list_size++] = n1;
                    }else if(fix_type[n1] == 0){
                        heap_fscore. add_node(n1);
                    }
                    score[n1]--;
                    freq_score[n1]-=freq[n1];
                    for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                        uint32_t n2 = adj[n1][n2_idx];
                        if(cds[n2]==0){
                            score[n2]--;
                            freq_score[n2]-=freq[n1];
                            if(fix_type[n2]!=2 && is_grey[n2])
                                heap_fscore.conduct_down(heap_fscore.idx_in_my_heap[n2]);
                        }
                    }
                }else if(dom_num[n1]==2){
                    uint32_t n1_dom_by = dom_by[n1];
                    ++score[n1_dom_by];
                    freq_score[n1_dom_by] += freq[n1];
                }
            }
            if(add_node != rand_pick_root){
                score[add_node] = -v_score;
                freq_score[add_node] = -v_fscore;
            }else{
                score[add_node] = -v_score-1;
                freq_score[add_node] = -v_fscore - freq[add_node];
            }
            
        }
    
    
    for(uint32_t vertex=1;vertex<=v_nums;++vertex){
        if(is_grey[vertex]==1){
            grey_vec[grey_vec_size] = vertex;
            idx_in_grey_vec[vertex] = grey_vec_size++;
        }else if(cds[vertex]==0){
            white_vec[white_vec_size] = vertex;
            idx_in_white_vec[vertex]  = white_vec_size++;
        }
    }
    undom_size = v_nums-dominated_sz;
    delete [] wait_add_list;
    // if(_opt_verbosity>0 && should_print_soln){
    //     cout<<"c HDC(Restart) with upperbound="<<upper_bound<<" ";
    //     cout<<"undom_size="<<undom_size<<" ";
    //     cout<<"restart time="<<++reconstruct_ct<<endl;
    // }
//    check_tmp_soln_with_score();
}


// a method used for test and comparson, 
uint32_t FastCDS::construct_by_outer(){
    cout << "c Do not use construct by outer! it only use for test at begining"<<endl;
    exit(0);
    /*
    ifstream fin(outer_cds_filename);
    uint32_t pick_node;
    set<uint32_t> node_set;
    while(fin>>pick_node){node_set.insert(pick_node);}
    fin.close();
    
    
    fix_points();
    fill_n(cds,         v_nums+1,   0);
    fill_n(is_grey,     v_nums+1,   0);
    uint32_t rand_root = *node_set.begin();
    is_grey[rand_root] = 1;
    for(uint32_t idx=1;idx<=v_nums;++idx) score[idx] = degree[idx]+1;
    score[rand_root]--;
    for(uint32_t n1_idx=0;n1_idx<degree[rand_root];++n1_idx) score[adj[rand_root][n1_idx]]--;
    
    
    int32_t     add_node            = -1;
    uint32_t    dominated_sz        = 1;
                cds_size            = 0;
                white_vec_size      = 0;
                grey_vec_size       = 0;
    uint32_t    wait_add_list_size  = 0;
    uint32_t    *wait_add_list      = new uint32_t[fixed_num];
    wait_add_list[wait_add_list_size++]=rand_root;
    
    while(dominated_sz < v_nums){
        add_node = wait_add_list[--wait_add_list_size];
        
        cds_vec[cds_size]   = add_node;
        idx_in_cds_vec[add_node]  = cds_size++;
        dominated_sz        += score[add_node];
        is_grey[add_node]   = 0;
        cds[add_node]       = 1;
        uint32_t v_score    = score[add_node];
        if(++dom_num[add_node]==2) score[dom_by[add_node]]++;
        else if(dom_num[add_node]==1) dom_by[add_node] = add_node;
        for(uint32_t n1_idx=0;n1_idx<degree[add_node];++n1_idx){
            uint32_t n1 = adj[add_node][n1_idx];
            dom_num[n1]++;
            if(dom_num[n1]==1){
                //white -> grey
                is_grey[n1] = 1;
                dom_by[n1]  = add_node;
                if(cds[n1]==0 && node_set.find(n1)!=node_set.end()){
                    wait_add_list[wait_add_list_size++] = n1;
                }
                score[n1]--;
                for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                    uint32_t n2 = adj[n1][n2_idx];
                    if(cds[n2]==0){
                        score[n2]--;
                    }
                }
            }else if(dom_num[n1]==2){
                score[dom_by[n1]]++;
            }
        }
        if(add_node != max_degree_vertex){
            score[add_node] = -v_score;
        }else{
            score[add_node] = -v_score-1;
        }
    }
    for(uint32_t vertex=1;vertex<=v_nums;++vertex){
        if(cds[vertex]==0){
            grey_vec[grey_vec_size] = vertex;
            idx_in_grey_vec[vertex] = grey_vec_size++;
        }
    }
    undom_size = 0;
    delete [] wait_add_list;
    copy(score, score+v_nums+1, freq_score);
    update_best_soln();
    if(_opt_verbosity>0){
        cout<<"c Initial cds_size = "<<best_cds_size<<endl;
        double runtime = getCpuTime();
        cout<<"c --TIME-- construct(Outer) use: "<<runtime-last_time_point<<" s"<<endl;
        last_time_point = runtime;
    }
    return cds_size;
    */
}


void FastCDS::add_to_cds(uint32_t vertex){
    pick_add_num[vertex]++;
    pick_swap_num[vertex]++;
    // cout<<" "<<vertex<<endl;

    cds[vertex]             =   1;
    cds_vec[cds_size]       =   vertex;
    idx_in_cds_vec[vertex]  =   cds_size++;
    // the vertex add to cds only can be the grey points.
    // grey -> black
    is_grey[vertex]             =   0;
    uint32_t idx_grey           =   idx_in_grey_vec[vertex];
    uint32_t last_node          =   grey_vec[--grey_vec_size];
    idx_in_grey_vec[last_node]  =   idx_grey;
    grey_vec[idx_grey]          =   last_node;
    int32_t vf_score            =   freq_score[vertex];
    int32_t v_score             =   score[vertex];
    undom_size                  -=  score[vertex];
    time_stamp[vertex]          = steps;
    
    if(++dom_num[vertex]==2){
        uint32_t v_dom_by = dom_by[vertex];
        ++score[v_dom_by];
        freq_score[v_dom_by] += freq[vertex];
    }else if(dom_num[vertex]==1) dom_by[vertex]=vertex;
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1 = adj[vertex][n1_idx];
        dom_num[n1]++;
        if(dom_num[n1] == 1){
            // white -> grey
            uint32_t  idx_white         = idx_in_white_vec[n1];
            uint32_t  last_white        = white_vec[--white_vec_size];
            idx_in_white_vec[last_white]= idx_white;
            white_vec [idx_white]       = last_white;
            is_grey[n1]                 = 1;
            grey_vec[grey_vec_size]     = n1;
            idx_in_grey_vec[n1]         = grey_vec_size++;
            dom_by[n1]      = vertex;
            score[n1]--;
            freq_score[n1] -= freq[n1];
            for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                uint32_t n2 = adj[n1][n2_idx];
                if(cds[n2] == 0){
                    score[n2]--;
                    freq_score[n2] -= freq[n1];
                }
            }
        }else if(dom_num[n1] == 2){
            // is grey node and the related black node will be effected.
            uint32_t n1_dom_by = dom_by[n1];
            ++score[n1_dom_by];
            freq_score[n1_dom_by] += freq[n1];
        }
    }
    score[vertex]       = - v_score;
    freq_score[vertex]  = - vf_score;
    
    //cc _ tabu
#ifdef ORIGIN_CC_TABU_METHOD
    uint8_t cc_tmp = conf_change[vertex];
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1 = adj[vertex][n1_idx];
        conf_change[n1] = 1;
        for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx)
            conf_change[adj[n1][n2_idx]] = 1;
    }
    conf_change[vertex] = cc_tmp;
    tabu_remove[vertex] = (uint64_t)(5+rand()%10);
#else
    uint8_t cc_tmp = conf_change[vertex];
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1 = adj[vertex][n1_idx];
        conf_change[n1] = 1;
        for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx)
            conf_change[adj[n1][n2_idx]] = 1;
    }
    conf_change[vertex] = cc_tmp;
    tabu_remove[vertex] = steps + 3;
#endif
}

void FastCDS::remove_from_cds(uint32_t vertex){
    pick_rm_num[vertex]++;
    pick_swap_num[vertex]++;
    // cout<<vertex<<" ";

    uint32_t v_idx          = idx_in_cds_vec[vertex];
    uint32_t last_v         = cds_vec[--cds_size];
    cds_vec[v_idx]          = last_v;
    idx_in_cds_vec[last_v]  = v_idx;
    cds[vertex]             = 0;
    time_stamp[vertex]      = steps;
    
    //from black -> grey
    is_grey[vertex]         =   1;
    grey_vec[grey_vec_size] =   vertex;
    idx_in_grey_vec[vertex] =   grey_vec_size++;
    int32_t v_score         =   score[vertex];
    int32_t vf_score        =   freq_score[vertex];
    undom_size              -=  score[vertex];      //black node is negative
    if(--dom_num[vertex]==1){
        //vertex can only be grey by the remove strategy
        uint32_t v_dom_by = -1;
        for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
            uint32_t n1=adj[vertex][n1_idx];
            if(cds[n1]==1){
                v_dom_by = n1;
                dom_by[vertex] = n1;
                break;
            }
        }
        --score[v_dom_by];
        freq_score[v_dom_by] -= freq[vertex];
    }
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1 = adj[vertex][n1_idx];
        dom_num[n1]--;
        if(dom_num[n1]==0){
            //from grey to white
            is_grey[n1] = 0;
            uint32_t    last_grey_point     = grey_vec[--grey_vec_size];
            uint32_t    idx_grey_point      = idx_in_grey_vec[n1];
            idx_in_grey_vec[last_grey_point]= idx_grey_point;
            grey_vec[idx_grey_point]        = last_grey_point;
            white_vec[white_vec_size]       = n1;
            idx_in_white_vec[n1]            = white_vec_size++;
            score[n1]++;
            freq_score[n1]+=freq[n1];
            for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                uint32_t n2 = adj[n1][n2_idx];
                if(cds[n2]==0){
                    score[n2]++;
                    freq_score[n2] += freq[n1];
                }
            }
            
        }else if(dom_num[n1]==1){
            uint32_t n1_dom_by = 0;
            if(cds[n1]==1){
                n1_dom_by = n1;
                dom_by[n1] = n1;
            }else{
                for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx){
                    uint32_t n2=adj[n1][n2_idx];
                    if(cds[n2]==1){
                        n1_dom_by = n2;
                        dom_by[n1] = n2;
                        break;
                    }
                }
            }
            --score[n1_dom_by];
            freq_score[n1_dom_by] -= freq[n1];
        }
    }
    score[vertex]       = - v_score;
    freq_score[vertex]  = - vf_score;
    
    //cc && tabu
#ifdef ORIGIN_CC_TABU_METHOD
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1 = adj[vertex][n1_idx];
        for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx)
            conf_change[adj[n1][n2_idx]] = 1;
        conf_change[n1] = 1;
    }
    conf_change[vertex] = 0;
#else
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1 = adj[vertex][n1_idx];
        for(uint32_t n2_idx=0;n2_idx<degree[n1];++n2_idx)
            conf_change[adj[n1][n2_idx]] = 1;
        conf_change[n1] = 1;
    }
    conf_change[vertex] = 0;
    tabu_add[vertex]    = steps + 1;
#endif
}

void FastCDS::update_leaf_add(uint32_t vertex){
    // try to link the leaf on a non-leaf tree node.
    uint32_t v_father           = 0;
    uint32_t max_tree_degree    = 0;
    pick_cand_size              = 0;
    for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
        uint32_t n1     = adj[vertex][n1_idx];
        uint32_t tmp_d  = son_num[n1];
//        if(cds[n1]==1 && tmp_d >= max_tree_degree){
//            max_tree_degree = tmp_d;
//            v_father        = n1;
//        }
        if(cds[n1]==1){
            if(tmp_d > max_tree_degree){
                max_tree_degree = tmp_d;
                pick_cand_size  = 1;
                pick_cand[0]    = n1;
            }else if(tmp_d == max_tree_degree){
                pick_cand[pick_cand_size++] = n1;
            }
        }
    }
    v_father = pick_cand[rand()%pick_cand_size];
    if(son_num[v_father]++ == 0){
        is_leaf[v_father]           = 0;
        uint32_t last_node          = leaf_vec[--leaf_vec_size];
        uint32_t idx_node           = idx_in_leaf_vec[v_father];
        idx_in_leaf_vec[last_node]  = idx_node;
        leaf_vec[idx_node]          = last_node;
    }
    father[vertex]              = v_father;
    son_num[vertex]             = 0;
    is_leaf[vertex]             = 1;
    leaf_vec[leaf_vec_size]     = vertex;
    idx_in_leaf_vec[vertex]     = leaf_vec_size++;
    
}

void FastCDS::update_leaf_remove(uint32_t vertex){
    is_leaf[vertex]             = 0;
    uint32_t last_node          = leaf_vec[--leaf_vec_size];
    uint32_t idx_node           = idx_in_leaf_vec[vertex];
    idx_in_leaf_vec[last_node]  = idx_node;
    leaf_vec[idx_node]          = last_node;
    
    uint32_t v_father = father[vertex];
    if(--son_num[v_father] == 0){
        is_leaf[v_father]           = 1;
        leaf_vec[leaf_vec_size]     = v_father;
        idx_in_leaf_vec[v_father]   = leaf_vec_size++;
    }
}


//For the main algorithm
//================================================
void FastCDS::qsort_by_score(uint32_t *a, uint32_t size){
    // the size is usually less than 100. so use insert sort instead.
    uint32_t    biggest_idx;
    int32_t     biggest_score;
    for(uint32_t ct=0;ct<size;++ct){
        biggest_idx = ct;
        biggest_score = freq_score[a[ct]];
        for(uint32_t inner=ct+1;inner<size;++inner){
            int32_t inner_score = freq_score[a[inner]];
            if(biggest_score < inner_score){
                biggest_score   = inner_score;
                biggest_idx     = inner;
            }
        }
        swap(a[ct],a[biggest_idx]);
    }
}


uint32_t FastCDS::pick_add(){ //too costy when degree is very big
    //add from grey nodes neighboring with at least one white nodes.
    uint32_t    add_node    = 0;
    int32_t     b_score     = INT32_MIN;
    pick_cand_size          = 0;
    // add a random white grey node
#ifdef RAND_BREAKING_TIE
    break_tie_vec_sz        = 0;
#endif
#ifdef ORIGIN_PICK_REMOVE
    for(uint32_t pick_idx=0;pick_idx<white_vec_size;++pick_idx){
        uint32_t pick_white = white_vec[pick_idx];
#else
    uint32_t pick_time  = min(white_vec_size,_opt_bms_pick_add);
    for(uint32_t pick_ct=0;pick_ct<pick_time;++pick_ct){
        uint32_t pick_white = white_vec[rand()%white_vec_size];
#endif
    
        uint32_t n1_sz      = degree[pick_white];
        for(uint32_t n1_idx=0;n1_idx<n1_sz;n1_idx++){
            uint32_t pick_node = adj[pick_white][n1_idx];
            if(is_grey[pick_node]){
                pick_cand[pick_cand_size++] = pick_node;
#ifdef REASONING_FORCE
                if(fix_type[pick_node] == 2) continue;
#endif
#ifdef REASONING_MORE_PROB
                if(fix_type[pick_node]==1) pick_cand[pick_cand_size++] = pick_node;
#endif
#ifdef ORIGIN_CC_TABU_METHOD
                if(conf_change[pick_node]==1){
#else
                if(conf_change[pick_node]==1 && steps > tabu_add[pick_node]){
#endif
                
                    int32_t cur_score  = freq_score[pick_node];
                    if(cur_score > b_score){
#ifdef RAND_BREAKING_TIE
                        break_tie_vec_sz = 1;
                        break_tie_vec[0] = pick_node;
#endif
                        add_node    = pick_node;
                        b_score     = cur_score;
                    }else if(cur_score == b_score){
#ifdef USE_SAFETY_AGE
                        float pick_node_score = dom_num[pick_node]*(steps-time_stamp[pick_node]);
                        float add_node_score  = dom_num[add_node]*(steps-time_stamp[add_node]);
                        if(pick_node_score > add_node_score){
                            add_node = pick_node;
                        }
#else
    #ifdef RAND_BREAKING_TIE
        #ifdef USE_SAFETY
                        if(dom_num[pick_node] > dom_num[add_node]){
                            break_tie_vec_sz = 1;
                            break_tie_vec[0] = pick_node;
                            add_node         = pick_node;
                        }else if(dom_num[pick_node] == dom_num[add_node]){
                            break_tie_vec[break_tie_vec_sz++] = pick_node;
                        }
        #else
                        break_tie_vec[break_tie_vec_sz++] = pick_node;
        #endif
    #else
        #ifdef USE_SAFETY
                        if((dom_num[pick_node] > dom_num[add_node]) ||
                           ((dom_num[pick_node] == dom_num[add_node])&&(time_stamp[pick_node] < time_stamp[add_node]))
                           ){
        #else
                        if(time_stamp[pick_node] < time_stamp[add_node]){
        #endif
                            add_node = pick_node;
                        }
    #endif
#endif
                    }
                }
            }
        }
    }
    if(add_node != 0){
#ifdef RAND_BREAKING_TIE
        return break_tie_vec[rand()%break_tie_vec_sz];
#else
        return add_node;
#endif
    }else if(pick_cand_size>0){
        return pick_cand[rand()%pick_cand_size];
    }else{
        return grey_vec[rand()%grey_vec_size];
    }

}

uint32_t FastCDS::pick_add_small(){
    uint32_t    add_node    = 0;
    int32_t     b_score     = INT32_MIN;
    pick_cand_size          = 0;
    // add a random white grey node
    uint32_t pick_time  = min((int32_t)white_vec_size,5);
    for(uint32_t pick_ct=0;pick_ct<pick_time;++pick_ct){
        uint32_t pick_white = white_vec[rand()%white_vec_size];
        uint32_t n1_sz      = degree[pick_white];
        for(uint32_t n1_idx=0;n1_idx<n1_sz;n1_idx++){
            uint32_t pick_node = adj[pick_white][n1_idx];
            if(is_grey[pick_node] && fix_type[pick_node]!=2){
                pick_cand[pick_cand_size++] = pick_node;
                if(conf_change[pick_node]==1 && steps > tabu_add[pick_node]){
                    int32_t cur_score  = freq_score[pick_node];
                    if(cur_score > b_score){
                        add_node    = pick_node;
                        b_score     = cur_score;
                    }else if(cur_score == b_score){
                        if((dom_num[pick_node] > dom_num[add_node]) ||
                           ((dom_num[pick_node] == dom_num[add_node])&&(time_stamp[pick_node] < time_stamp[add_node]))
                           ){
                            add_node = pick_node;
                        }
                    }
                }
            }
        }
    }
    if(add_node != 0){
        return add_node;
    }else if(pick_cand_size>0){
        return pick_cand[rand()%pick_cand_size];
    }else{
        return grey_vec[rand()%grey_vec_size];
    }
}


uint32_t FastCDS::pick_remove(){
    // pick node from leaf_node.
    // all remove nodes are negative,
    // find the max score nodes.
    // all leaf nodes can not be cut points, so it can br removed unless be CCed or TABUed.
#ifdef RAND_BREAKING_TIE
    break_tie_vec_sz=0;
#endif
    uint32_t    remove_node     = 0;
    int32_t     b_freq_score    = INT32_MIN;
    pick_cand_size       = 0;
    for(uint32_t ct=0;ct<_opt_bms_pick_remove;++ct){
        uint32_t pick_node;
        if(TBC){
            pick_node = leaf_vec[rand()%leaf_vec_size];
        }else{
            pick_node = cds_vec[rand()%cds_size];
            if(is_cut[pick_node]==1) continue;
        }
        pick_cand[pick_cand_size++] = pick_node;
#ifdef REASONING_FORCE
        if(fix_type[pick_node] == 1) continue;
#endif 
#ifdef REASONING_MORE_PROB
        if(fix_type[pick_node] == 2) pick_cand[pick_cand_size++] = pick_node;
#endif        
        if(steps > tabu_remove[pick_node]){
            int32_t cur_score  = freq_score[pick_node];
            if(cur_score > b_freq_score){
#ifdef RAND_BREAKING_TIE
                break_tie_vec_sz = 1;
                break_tie_vec[0] = pick_node;              
#endif          
                remove_node     = pick_node;
                b_freq_score    = cur_score;
            }else if(cur_score == b_freq_score){
#ifdef USE_SAFETY_AGE
                float pick_node_subscore   = dom_num[pick_node]*(time_stamp[pick_node]+1);
                float remove_node_subscore = dom_num[remove_node]*(time_stamp[remove_node]+1);
                // cout<<pick_node_subscore<<" "<<remove_node_subscore<<endl;
                if(pick_node_subscore<remove_node_subscore){
                    remove_node  = pick_node;
                }
#else
    #ifdef RAND_BREAKING_TIE
        #ifdef USE_SAFETY
                if( dom_num[pick_node]<dom_num[remove_node] ){
                    remove_node = pick_node;
                    break_tie_vec_sz = 1;
                    break_tie_vec[0] = pick_node;     
                }else if (dom_num[pick_node] == dom_num[remove_node]){
                    break_tie_vec[break_tie_vec_sz++] = pick_node;
                }
        #else
                break_tie_vec[break_tie_vec_sz++] = pick_node;
        #endif
    #else
        #ifdef USE_SAFETY
                if( (dom_num[pick_node]<dom_num[remove_node])
                   ||(dom_num[pick_node]==dom_num[remove_node] && time_stamp[pick_node]<time_stamp[remove_node])
                ){
        #else
                if( time_stamp[pick_node]<time_stamp[remove_node] ){
        #endif
                    remove_node  = pick_node;
                }
    #endif
#endif
            }
        }
    }
    if(remove_node!=0){
#ifdef RAND_BREAKING_TIE
        return break_tie_vec[rand()%break_tie_vec_sz];
#else
        return remove_node;
#endif
    }else if(pick_cand_size>0){
#ifdef ORIGIN_PICK_REMOVE
        qsort_by_score(pick_cand,pick_cand_size);
        uint32_t pick_cand_mod = pick_cand_size*0.4+1;
        return pick_cand[rand()%pick_cand_mod];
#else     
        return pick_cand[rand()%pick_cand_size];
#endif
        // for(uint32_t node_idx=0;node_idx<pick_cand_size;++node_idx){
        //     uint32_t pick_node = pick_cand[node_idx];
        //     int32_t cur_score  = freq_score[pick_node];
        //     if(cur_score > b_freq_score){
        //         remove_node     = pick_node;
        //         b_freq_score    = cur_score;
        //     }else if(cur_score == b_freq_score){
        //         if( (dom_num[pick_node]<dom_num[remove_node])
        //            ||(dom_num[pick_node]==dom_num[remove_node] && time_stamp[pick_node]<time_stamp[remove_node])
        //         ){
        //             remove_node  = pick_node;
        //         }
        //     }
        // }
        // return remove_node;
    }else{
        if(TBC){
            return leaf_vec[rand()%leaf_vec_size];
        }else{
            while(true){
                remove_node = cds_vec[rand()%cds_size];
                if(is_cut[remove_node]==0) return remove_node;
            }
        }
    }
}


uint32_t FastCDS::pick_remove_small(){
    // pick node from leaf_node.
    // all remove nodes are negative,
    // find the max score nodes.
    // all leaf nodes can not be cut points, so it can br removed unless be CCed or TABUed.
    uint32_t    remove_node     = 0;
    int32_t     b_freq_score    = INT32_MIN;
    pick_cand_size              = 0;
    
    uint32_t pick_idx   = 0;
    uint32_t pick_bound;
    if(TBC)
        pick_bound      =   leaf_vec_size;
    else
        pick_bound      =   cds_size;
    
    for(;pick_idx<pick_bound;pick_idx++){
        uint32_t pick_node;
        if(TBC){
            pick_node = leaf_vec[pick_idx];
        }else{
            pick_node = cds_vec[pick_idx];
            if(is_cut[pick_node]==1) continue;
        }
        pick_cand[pick_cand_size++] = pick_node;
        if(steps > tabu_remove[pick_node]){
            int32_t cur_score  = freq_score[pick_node];
            if(cur_score > b_freq_score){
                remove_node     = pick_node;
                b_freq_score    = cur_score;
            }else if(cur_score == b_freq_score){
                if( (dom_num[pick_node]<dom_num[remove_node])
                   ||(dom_num[pick_node]==dom_num[remove_node] && time_stamp[pick_node]<time_stamp[remove_node])
                ){
                    remove_node  = pick_node;
                }
            }
        }

    }
    if(remove_node!=0){
        return remove_node;
    }
    else{
        // cout<<"stucking remove"<<endl;
        return pick_cand[rand()%pick_cand_size];
    }
    
    
}



void FastCDS::update_freq_scores(){
    for(uint32_t idx=0;idx<white_vec_size;idx++){
        uint32_t white_node = white_vec[idx];
        freq[white_node]++;
        freq_score[white_node]++;
        uint32_t white_degree = degree[white_node];
        for(uint32_t n1_idx=0;n1_idx<white_degree;++n1_idx){
            uint32_t n1 = adj[white_node][n1_idx];
            freq_score[n1]++;
        }
    }
}

void FastCDS::smooth_weight(){
    smooth_weight_ct++;
    for(uint32_t vertex=1;vertex<=v_nums;++vertex){
        uint32_t new_freq = 0.3*freq[vertex];
        if(new_freq<1) new_freq = 1;
        uint32_t loose = freq[vertex]-new_freq;
        if(loose>0){
            freq[vertex]    =   new_freq;
            if(cds[vertex]==0 && is_grey[vertex]==0){
                //white node
                freq_score[vertex] -= loose;
                for(uint32_t n1_idx=0;n1_idx<degree[vertex];++n1_idx){
                    freq_score[adj[vertex][n1_idx]] -= loose;
                }
            }else if(dom_num[vertex]==1){
                //grey node only domed by one black node;
                freq_score[dom_by[vertex]] += loose;
            }
            total_weight    -=  loose;
        }
    }
}

uint32_t FastCDS::do_LS_SUB(){
    TBC = false;
    uint32_t    improve_this_turn   = 0;
    uint32_t    no_improve_steps    = 0;
    uint32_t    step_this_turn      = 1;
    //preprocess



    while(true){
        
        if(cds_size<5000 && increase_SUB_flag){
            increase_SUB_flag           =   false;
            min_base_SUB                *=  10;
            max_steps_SUB               *=  10;
            tmp_max_steps_noimpr_SUB    *=  10;
        }
        
        if(step_this_turn % try_step_SUB==0){
            double tmp_runtime = getCpuTime();
            if(should_print_soln==false && tmp_runtime>(_opt_cutoff_time-_opt_print_soln_gap)){
                cout<<"o "<<best_cds_size<<" "<<best_cds_time<<endl;
                should_print_soln = true;
            }
            if(tmp_runtime>_opt_cutoff_time)
                return 0;
            if(_opt_for_small_ins==1){
                if(step_this_turn > max_steps_SUB){
                    if(improve_this_turn>10){
                        tmp_max_steps_noimpr_SUB -= min_base_SUB;
                        if(tmp_max_steps_noimpr_SUB < min_base_SUB){
                            tmp_max_steps_noimpr_SUB = min_base_SUB;
                        }
                    }
                    return best_cds_size;
                }
                if(no_improve_steps>tmp_max_steps_noimpr_SUB){
                    tmp_max_steps_noimpr_SUB += min_base_SUB;
                    if(tmp_max_steps_noimpr_SUB > max_steps_SUB){
                        tmp_max_steps_noimpr_SUB = max_steps_SUB;
                    }
                    return best_cds_size;
                }
            }else{
                if (step_this_turn > max_steps_SUB 
                || no_improve_steps>tmp_max_steps_noimpr_SUB){
                    return best_cds_size;
                }
            }
            
        }
        if(undom_size == 0){
            update_best_soln();
            if(cds_size == 2) return 2;
            renew_cut();
            uint32_t rm_node;
            if(_opt_for_small_ins==1)
                rm_node = pick_remove_small();
            else
                rm_node = pick_remove();
#ifdef COUNT_PICK_REMOVE
        SUB_pick_num += 1;
        ct_remove_SUB[rm_node] += 1;
#endif
#ifdef ANALYZE_HDC
        if(SUB_remove.find(rm_node) == SUB_remove.end()){
            SUB_remove[rm_node]   = 1;
        }else SUB_remove[rm_node] +=1;
#endif
            // cout<<"rm:"<<rm_node<<endl;
            remove_from_cds(rm_node);
            conf_change[rm_node]=   1;
            no_improve_steps    =   0;
            improve_this_turn   ++;
            continue;
        }
        

        //remove
        renew_cut();
        uint32_t rm_node;
        if(_opt_for_small_ins==1)
            rm_node = pick_remove_small();
        else
            rm_node = pick_remove();
        remove_from_cds(rm_node);

        //add_new
        uint32_t add_node;
        if(_opt_for_small_ins==1)
            add_node = pick_add_small();
        else
            add_node = pick_add();
        add_to_cds(add_node);

#ifdef COUNT_PICK_REMOVE
        SUB_pick_num += 1;
        ct_remove_SUB[rm_node] += 1;
#endif

#ifdef ANALYZE_HDC
        if(SUB_add.find(add_node) == SUB_add.end()){
            SUB_add[add_node]   = 1;
        }else SUB_add[add_node] +=1;
        if(SUB_remove.find(rm_node) == SUB_remove.end()){
            SUB_remove[rm_node]   = 1;
        }else SUB_remove[rm_node] +=1;
#endif
        

        
        update_freq_scores();
        
        //smooth weight
        total_weight += undom_size;
        if(total_weight > weight_threshold){
            smooth_weight();
        }
        
        //update score of undomed variables;
        
        steps++;
        SUB_steps++;
        no_improve_steps++;
        step_this_turn++;
    }
    



    return best_cds_size;
}


uint32_t FastCDS::do_LS_TBC_without_HDC(){
    TBC = true;
    uint32_t    improve_this_turn   = 0;
    uint32_t    no_improve_steps    = 0;
    uint32_t    step_this_turn      = 1;
    //preprocess
    construct_tree_random_based();
       
    while(true){
        if(step_this_turn % try_step_TBC == 0){
            double tmp_runtime = getCpuTime();
            if(should_print_soln==false && tmp_runtime>(_opt_cutoff_time-_opt_print_soln_gap)){
                cout<<"o "<<best_cds_size<<" "<<best_cds_time<<endl;
                should_print_soln = true;
            }
            if(tmp_runtime>_opt_cutoff_time)
                return 0;
        }
        if(undom_size == 0){
            update_best_soln();
            if(cds_size == 2) return 2;
            uint32_t rm_node;
            rm_node = pick_remove();
            
            remove_from_cds(rm_node);
            update_leaf_remove(rm_node);
            conf_change[rm_node]    = 1;
            no_improve_steps        = 0;
            improve_this_turn++;
            continue;
        }

        //remove
        uint32_t rm_node;
        rm_node = pick_remove();
        remove_from_cds(rm_node);
        update_leaf_remove(rm_node);

        //add
        uint32_t add_node;
        add_node = pick_add();
        add_to_cds(add_node);
        update_leaf_add(add_node);
        
        update_freq_scores();
        
        //smooth weight
        total_weight += undom_size;
        if(total_weight > weight_threshold){
            smooth_weight();
        }
        //update score of undomed variables;
        steps++;
        no_improve_steps++;
        step_this_turn++;
        if(no_improve_steps % _opt_reconstruct_gap==0){
            construct_tree_random_based();
        }
                
    }
    
    
    return best_cds_size;
}


uint32_t FastCDS::do_LS_TBC(){
    TBC = true;
    uint32_t    improve_this_turn   = 0;
    uint32_t    no_improve_steps    = 0;
    uint32_t    step_this_turn      = 1;
    //preprocess
    
    if(_opt_for_small_ins==1){
        construct_tree_random_based();
    }else{
        construct_tree_random_based();
        // construct_tree_random_based_for_big();
        // construct_tree_hueristic_based();
        // construct_tree_subgraph_based();
    }
    while(true){
        if(step_this_turn % try_step_TBC == 0){
            double tmp_runtime = getCpuTime();
            if(should_print_soln==false && tmp_runtime>(_opt_cutoff_time-_opt_print_soln_gap)){
                cout<<"o "<<best_cds_size<<" "<<best_cds_time<<endl;
                should_print_soln = true;
            }
            if(tmp_runtime>_opt_cutoff_time)
                return 0;
            if(step_this_turn > max_steps_TBC){
                if(improve_this_turn>10){
                    tmp_max_steps_noimpr_TBC -= min_base_TBC;
                    if(tmp_max_steps_noimpr_TBC < min_base_TBC){
                        tmp_max_steps_noimpr_TBC = min_base_TBC;
                    }
                }
                return best_cds_size;
            }
            if(no_improve_steps > tmp_max_steps_noimpr_TBC){
                tmp_max_steps_noimpr_TBC += min_base_TBC;
                if(tmp_max_steps_noimpr_TBC > max_steps_TBC){
                    tmp_max_steps_noimpr_TBC = max_steps_TBC;
                }
                return best_cds_size;
            }
        }
        if(undom_size == 0){
            update_best_soln();
            if(cds_size == 2) return 2;
            uint32_t rm_node;
            if(_opt_for_small_ins==1)
                rm_node = pick_remove_small();
            else
                rm_node = pick_remove();
#ifdef COUNT_PICK_REMOVE
            TBC_pick_num += 1;
            ct_remove_TBC[rm_node] += 1;
#endif
#ifdef ANALYZE_HDC
            if(TBC_remove.find(rm_node) == TBC_remove.end()){
                TBC_remove[rm_node]   = 1;
            }else TBC_remove[rm_node] +=1;
#endif
            remove_from_cds(rm_node);
            update_leaf_remove(rm_node);
            conf_change[rm_node]    = 1;
            no_improve_steps        = 0;
            improve_this_turn++;
            continue;
        }



        //remove
        uint32_t rm_node;
        if(_opt_for_small_ins==1)
            rm_node = pick_remove_small();
        else
            rm_node = pick_remove();
        remove_from_cds(rm_node);
        update_leaf_remove(rm_node);
        


        //add
        uint32_t add_node;
        if(_opt_for_small_ins==1)
            add_node = pick_add_small();
        else
            add_node = pick_add();
        add_to_cds(add_node);
        update_leaf_add(add_node);

#ifdef COUNT_PICK_REMOVE
        TBC_pick_num += 1;
        ct_remove_TBC[rm_node] += 1;
#endif
#ifdef ANALYZE_HDC
        if(TBC_add.find(add_node) == TBC_add.end()){
            TBC_add[add_node]   = 1;
        }else TBC_add[add_node] +=1;
        if(TBC_remove.find(rm_node) == TBC_remove.end()){
            TBC_remove[rm_node]   = 1;
        }else TBC_remove[rm_node] +=1;
#endif
        
        
        
        update_freq_scores();
        
        //smooth weight
        total_weight += undom_size;
        if(total_weight > weight_threshold){
            smooth_weight();
        }
        //update score of undomed variables;
        steps++;
        no_improve_steps++;
        step_this_turn++;
#ifndef ANALYZE_HDC
        if(no_improve_steps % _opt_reconstruct_gap==0){
            reconstruct_num++;
            if(_opt_for_small_ins==1){
                construct_tree_random_based();
            }else{
                construct_tree_random_based();
                // construct_tree_random_based_for_big();
                // construct_tree_hueristic_based();
                // construct_tree_subgraph_based();
            }
        }
#endif
    }
    
    
    return best_cds_size;
}


bool FastCDS::exec_only_TBC_plus(){
    uint32_t res;
    while(true){
        
        restarts++;
        res = do_LS_TBC_without_HDC();
        reconstruct_with_upper_bound(best_cds_size-1);
        if(res == 0) return false;
        if(res == 2){
            cout<<"c best soln found!"<<endl;
            break;
        }
    }
}

bool FastCDS::exec_for_classical(){
    uint32_t res;
    while(true){
        
        restarts++;
        double in_time = getCpuTime();
        res = do_LS_SUB();
        double SUB_used_time = getCpuTime()-in_time;
        SUB_sum_time += SUB_used_time;
        reconstruct_with_upper_bound(best_cds_size-1);
        if(res == 0) return false;
        if(res == 2){
            cout<<"c best soln found!"<<endl;
            break;
        }
    }
}

bool FastCDS::exec_for_massive(){
    uint32_t res;
    while(true){
        uint32_t before_cds = best_cds_size;
        if(_opt_only_SUB == 0){
            if(restarts++%2==0){
                double in_time = getCpuTime();
                res = do_LS_TBC();
                double TBC_used_time = getCpuTime()-in_time;
                TBC_sum_time += TBC_used_time;
            }else{
#ifdef ANALYZE_HDC 
                renew_cut();
                for(uint32_t idx=0;idx<cds_size;++idx){
                    uint32_t node = cds_vec[idx];
                    if(is_cut[node]==0){
                        before_SUB.insert(node);
                    }
                }
#endif          
                double in_time = getCpuTime();
                res = do_LS_SUB();
                double SUB_used_time = getCpuTime()-in_time;
                SUB_sum_time += SUB_used_time;
#ifdef ANALYZE_HDC 
                renew_cut();
                for(uint32_t idx=0;idx<cds_size;++idx){
                    uint32_t node = cds_vec[idx];
                    if(is_cut[node]==0){
                        after_SUB.insert(node);
                    }
                }
#endif
            }
#ifdef ANALYZE_HDC
            uint32_t ct = restarts/2;
            if(restarts%2==0){
                if(ct > 0){
                    cout<<"debug_info "
                        // <<TBC_add.size()<<" "
                        <<"TBC_last  : "<<TBC_remove_last.size()<<" "
                        <<"TBC_remove: "<<TBC_remove.size()<<" "
                        <<"SUB_remove: "<<SUB_remove.size()<<" "
                        <<"before_SUB: "<<before_SUB.size()<<" "
                        <<"after_SUB : "<<after_SUB.size()<<" "
                        <<endl;
                    set<uint32_t> s_joint,s_unoin;
                    vector<pair<uint32_t,uint32_t>> sort_TBC;
                    vector<pair<uint32_t,uint32_t>> sort_SUB;
                    vector<uint32_t> SUB_position;
                    uint32_t smaller_sz;
                    // compare add
                    // s_joint.clear(),s_unoin.clear();
                    // sort_TBC.clear();
                    // for(auto item:TBC_add) sort_TBC.push_back(item);
                    // sort(sort_TBC.begin(),sort_TBC.end(),
                    //     [](pair<uint32_t,uint32_t> a, pair<uint32_t,uint32_t> b)->bool{
                    //     return (a.second > b.second);
                    // });
                    // // for(auto item : SUB_remove){cout<<"e "<<item.first<<" "<<item.second<<endl;}
                    // smaller_sz = min(SUB_add.size(),TBC_add.size())/2;
                    // for(auto item:SUB_add) s_unoin.insert(item.first);
                    // for(uint32_t idx=0;idx<smaller_sz;++idx){
                    //     uint32_t elem = sort_TBC[idx].first;
                    //     s_unoin.insert(elem);
                    //     if(SUB_add.find(elem)!=SUB_add.end())
                    //         s_joint.insert(elem);
                    // }
                    // cout<<"debug_add "
                    //     <<s_joint.size()<<" "
                    //     <<s_unoin.size()<<" "
                    //     <<(s_joint.size()+0.0)/(s_unoin.size())
                    //     <<endl;

                    // compare remove
                    s_joint.clear(),s_unoin.clear();
                    sort_TBC.clear(),sort_SUB.clear();
                    for(auto item:TBC_remove) sort_TBC.push_back(item);
                    sort(sort_TBC.begin(),sort_TBC.end(),
                        [](pair<uint32_t,uint32_t> a, pair<uint32_t,uint32_t> b)->bool{
                        return (a.second > b.second);
                    });
                    for(auto item:SUB_remove) sort_SUB.push_back(item);
                    sort(sort_SUB.begin(),sort_SUB.end(),
                        [](pair<uint32_t,uint32_t> a, pair<uint32_t,uint32_t> b)->bool{
                        return (a.second > b.second);
                    });
                    // for(auto item : sort_SUB){cout<<"e "<<item.first<<" "<<item.second<<endl;}
                    smaller_sz = min(SUB_remove.size(),TBC_remove.size());
                    SUB_position.clear();
                    for(auto item:sort_SUB){
                        uint32_t pos = 0;
                        for(uint32_t idx=0;idx<sort_TBC.size();++idx){
                            if(item.first == sort_TBC[idx].first){
                                pos = idx;
                                break;
                            }
                        }
                        SUB_position.push_back(pos);
                    }
                    cout<<"debug e ";
                    for(auto item:SUB_position) cout<<item<<" ";cout<<endl;

                    for(auto item:SUB_remove) s_unoin.insert(item.first);
                    for(uint32_t idx=0;idx<smaller_sz;++idx){
                        uint32_t elem = sort_TBC[idx].first;
                        s_unoin.insert(elem);
                        if(SUB_remove.find(elem)!=SUB_remove.end())
                            s_joint.insert(elem);
                    }
                    cout<<"debug_remove "
                        <<s_joint.size()<<" "
                        <<s_unoin.size()<<" "
                        <<(s_joint.size()+0.0)/(s_unoin.size())
                        <<endl;

                    
                    // compare last
                    s_joint.clear(),s_unoin.clear();
                    sort_TBC.clear();
                    for(auto item:TBC_remove_last) s_unoin.insert(item.first);
                    for(auto item:TBC_remove){
                        uint32_t elem = item.first;
                        s_unoin.insert(elem);
                        if(TBC_remove_last.find(elem)!=TBC_remove_last.end())
                            s_joint.insert(elem);
                    }
                    cout<<"debug_2TBC "
                        <<s_joint.size()<<" "
                        <<s_unoin.size()<<" "
                        <<(s_joint.size()+0.0)/(s_unoin.size())
                        <<endl;

                    // compare before after
                    s_joint.clear(),s_unoin.clear();
                    sort_TBC.clear();
                    for(auto item:before_SUB) s_unoin.insert(item);
                    for(auto elem:after_SUB){
                        s_unoin.insert(elem);
                        if(before_SUB.find(elem)!=before_SUB.end())
                            s_joint.insert(elem);
                    }
                    cout<<"debug_before_after "
                        <<s_joint.size()<<" "
                        <<s_unoin.size()<<" "
                        <<(s_joint.size()+0.0)/(s_unoin.size())
                        <<endl;



                    TBC_remove_last = TBC_remove;
                    TBC_add.clear();
                    TBC_remove.clear();
                    SUB_add.clear();
                    SUB_remove.clear();
                    
                }
            }

#endif
        }else{
            restarts++;
            res = do_LS_SUB();
        }

        if(res == 0) break;
        if(res == 2){
            cout<<"c best soln found!"<<endl;
            break;
        }
    }
}


//for the algorithm
bool FastCDS::exec(string filename){
    // reconstruct_signal = false;
    // ineffective_search_num = 0;
    start_time = chrono::steady_clock::now();
    SUB_sum_time = TBC_sum_time = 0;
    std::ios::sync_with_stdio(false);
    if(!build(filename)){
        cout<<"c --BUG-- build wrong"<<endl;
        return 10;
    }
    //construct
    if(outer_cds_filename==string("null")){
        if(_opt_bms_construct == 0)
            construct_1hop();
        else{
            construct_1hop_bms();
        }
    }else{
        construct_by_outer();
    }
    if(max_degree>=v_nums-1){
        cout<<"s 1 0.00";
        cout<<"c best soln found!"<<endl;
        cout<<"c end with construction"<<endl;
        return true;
    }

    if(_opt_only_init == 1){
        cout<<"c Cutoff immediately after initialize"<<endl;
        return false;
    }

    // exec_only_TBC_plus();
    if(_opt_for_small_ins == 1){
        exec_for_classical();
    }else{
        set_subgraph_data_structure();
        exec_for_massive();
    }

#ifdef COUNT_PICK_REMOVE
    cout<<"debug-info "<<SUB_pick_num<<" "<<TBC_pick_num<<endl;
    cout<<"debug-SUB-ct ";
    for(uint32_t var=1;var<=v_nums;++var) {
        cout<<ct_remove_SUB[var]<<" ";
        // if(ct_remove_SUB[var]>0) cout<<var<<","<<ct_remove_SUB[var]<<" ";
    }cout<<endl;
    cout<<"debug-TBC-ct ";
    for(uint32_t var=1;var<=v_nums;++var) {
        cout<<ct_remove_TBC[var]<<" ";
        // if(ct_remove_TBC[var]>0) cout<<var<<","<<ct_remove_TBC[var]<<" ";
    }cout<<endl;
#endif

    cout<<"c End with steps="<<steps<<"(SUB="<<SUB_steps
        <<"),restarts="<<restarts
        <<",smooth="<<smooth_weight_ct
        // <<",HDC(Reconstruct)="<<reconstruct_ct
        <<",TBC re-construction="<<reconstruct_num
        <<endl;
    cout<<"c SUB_time = "<<SUB_sum_time
        <<"  TBC_time = "<<TBC_sum_time
        <<endl;
    return true;
}

void show_basic_info(){
    cout<<"c ";
    for(uint32_t i=0;i<30;++i)cout<<"-";
    cout<<"  FastCDS  ";
    for(uint32_t i=0;i<30;++i)cout<<"-";cout<<endl;
    cout.setf(ios::left);
    cout<<setw(73)<< "c| FastCDS is based on NuCDS.";cout<<"|"<<endl;
    cout<<setw(73)<< "c| NuCDS: An Efficient Local Search Algorithm for Minimum Connected";cout<<"|"<<endl;
    cout<<setw(73)<< "c|        Dominating Set. IJCAI 2020.";cout<<"|"<<endl;
    cout<<setw(73)<< "c| CopyRight (c) 2020, Xindi Zhang, Bohan Li, Shaowei Cai, Yiyuan Wang.";cout<<"|"<<endl;
    cout<<setw(73)<< "c| Develop by Xindi Zhang (2020-06-05).";cout<<"|"<<endl;
    cout<<setw(73)<< "c| Email : dezhangxd@163.com.";cout<<"|"<<endl;
    cout<<setw(73)<< "c| This solver is under MIT License.";cout<<"|"<<endl;
    
    
    cout<<"c ";
    for(uint32_t i=0;i<71;++i)cout<<"-";cout<<endl;
}

int main(int argc, char *argv[]) {
    show_basic_info();
    parse_options(argc, argv);
    FastCDS fastCDS;
    string filename(argv[1]);
    fastCDS.exec(filename);
    fastCDS.print_soln(true);

    return 0;
}
