#ifndef VALOR_QC_GRAPH
#define VALOR_QC_GRAPH
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <algorithm>
#include <limits>


template <class V>
class Vertex{
public:
    Vertex(){}
	std::unordered_set<V *> edges;
	short tabu = 0;
	int dv = 0;
	bool covered = false;
	bool inactive = false;
};

template <class V>
class qcgraph{
	double lambda;
	double gamma;
    unsigned rest_tabu;
	std::unordered_map<V *,Vertex<V>> map;
    bool pairs_sorted = false;
    std::vector<std::pair<V * const,Vertex<V>> *> pair_ptrs; 
public:
	qcgraph(double lambda, double gamma, unsigned tabu);
    void add_vertex(V *sv);
	void add_edge(V &s1, V &s2);
	std::vector<V *> pop_quasi_clique(std::vector<V *> &component);
    void remove_all(std::vector<V *> vec);

    std::vector<std::vector<V *>> find_clusters();
    void dfs_step(std::vector<V *> &component, V *sv);
};


template<class V>
qcgraph<V>::qcgraph(double lambda, double gamma, unsigned tabu) : lambda(lambda), gamma(gamma), rest_tabu(tabu){}

template<class V>
void qcgraph<V>::add_vertex(V *sv){
    map[sv] = Vertex<V>();
}

template<class V>
void qcgraph<V>::add_edge(V &s1, V &s2){
    map[&s1].edges.insert(&s2);
    map[&s2].edges.insert(&s1);
}

template<class V>
void qcgraph<V>::dfs_step(std::vector< V*> &component, V *sv){
    map[sv].covered = true;
    component.push_back(sv);
    std::unordered_set<V *> edges = map[sv].edges;
    for( auto iter = edges.begin(); iter != edges.end(); iter++){
        
        if( !map[*iter].covered){
            dfs_step(component,*iter);
        }
    }
}

template<class V>
std::vector<std::vector<V *>> qcgraph<V>::find_clusters(){

    std::vector<std::vector<V *>> components;

    for(auto iter = map.begin(); iter != map.end(); iter++){
        if( !iter->second.covered){
            std::vector< V*> component;
            dfs_step(component,iter->first);
            components.push_back(component);
        }
    }

    //std::vector<std::vector<V *>> clusters;
    //for(auto iter = components.begin(); iter != components.end(); iter++){
    //    clusters.push_back(pop_quasi_clique(*iter));
    //}
    return components;
}

template<class V>
std::vector<V *> qcgraph<V>::pop_quasi_clique(std::vector< V*> &component){
    int e_prime = 0;
    int v_prime = 0;
    std::unordered_set<V *> clique_items;
    if(pair_ptrs.size() == 0){
        for(auto it = map.begin();it!=map.end();it++){
            it->second.dv = 1;
            pair_ptrs.push_back(&*it);
        }
    }

    if(!pairs_sorted){
        std::sort(pair_ptrs.begin(), pair_ptrs.end(), 
                [](const std::pair<V * const,Vertex<V>> *p1, const std::pair<V * const,Vertex<V>> *p2){
                return p1->second.edges.size() > p2->second.edges.size();
                });
    }

    bool added = true;
    bool removed = true;
    while( added || removed){
        added = false;
        removed = false;
        auto best_ptr = pair_ptrs.end() + 1;
        double max_score = std::numeric_limits<double>::min();
        for(auto pair_ptr = pair_ptrs.begin(); pair_ptr != pair_ptrs.end(); pair_ptr++){
            auto vertex = (*pair_ptr)->second;
            if( vertex.inactive){
                continue;
            }
            if( vertex.tabu > 0){
                vertex.tabu--;
            }
            if( vertex.tabu == 0 &&
                    vertex.dv > this->gamma * (v_prime -1) &&
                    (e_prime + vertex.dv >= (this->lambda * (v_prime * (v_prime-1)))/2.0)){
                double score = (2.0 * ( e_prime + vertex.dv)) / (v_prime * (v_prime +1));
                if( score > max_score){
                    max_score = score;
                    best_ptr = pair_ptr;
                }
            }        
        }
        if ( best_ptr != pair_ptrs.end() +1){

            
            clique_items.insert((*best_ptr)->first);
            v_prime++;

            added = true;
            auto &vertex = (*best_ptr)->second;
            vertex.inactive = true;
            for(auto sv : vertex.edges){
                if(clique_items.find(sv) != clique_items.end()){
                    e_prime++;
                }
                map[sv].dv++;
            }
            vertex.tabu = rest_tabu;

        }
        /*
        //Remove
        auto best_sv = clique_items.end();
        max_score = std::numeric_limits<double>::min();
        for( auto sv = clique_items.begin(); sv!=clique_items.end();sv++){
            auto vertex = map[*sv];
            if(vertex.tabu > 1){
                vertex.tabu--;

            }else if( vertex.dv < gamma * (v_prime -1) &&
                   (e_prime - vertex.dv > (lambda * ((v_prime - 2) * (v_prime -1))/2.0))){
                double score = (2.0 * (e_prime - vertex.dv)) / ((v_prime-1) *(v_prime-2));
                if( score > max_score){
                    max_score = score;
                    best_sv = sv;
                }
            } 
        }
        if( best_sv != clique_items.end()){
            removed = true;
            auto vertex = map[*best_sv];
            vertex.tabu = rest_tabu;
            vertex.inactive = false;
            for(auto sv:vertex.edges){
                if(clique_items.find(sv) != clique_items.end()){
                    e_prime--;
                }
                map[sv].dv--;
            }
            clique_items.erase( *best_sv);
            v_prime--;
            
        }
        */
    }
    std::vector<V *> cl_vector(clique_items.begin(),clique_items.end());
    remove_all(cl_vector);
    return cl_vector;
}


template<class V>
void qcgraph<V>::remove_all(std::vector<V *> vec){

    for( auto iter = vec.begin(); iter!= vec.end(); iter++){
        auto vertex = map[*iter];
        for( auto edger = vertex.edges.begin(); edger != vertex.edges.end();edger++){
            auto neighbour = map[*edger];
            neighbour.edges.erase(*iter);
        }
        map.erase(*iter);
    }
}
#endif
