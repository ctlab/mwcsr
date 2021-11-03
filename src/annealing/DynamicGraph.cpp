#include "include/DynamicGraph.h"

#include <cmath>
#include <utility>
#include <limits>
#include <iostream>

namespace dgraph {

    DynamicGraph::DynamicGraph(unsigned n) : n(n) {
        size = std::lround(std::ceil(std::log2(n)) + 1);
        for (unsigned i = 0; i < size; i++) {
            forests.emplace_back(n);
            adjLists.emplace_back();
            for (unsigned j = 0; j < n; j++){
                adjLists[i].push_back(new List());
            }
        }
    }

    DynamicGraph::~DynamicGraph() {
        for (unsigned i = 0; i < size; i++) {
            for (unsigned j = 0; j < n; j++) {
                ListIterator it = adjLists[i][j]->iterator();
                while (it.hasNext()) {
                    List* list = *it;
                    it++;
                    delete list->e();
                }
                delete *it;
            }
        }
    }

    EdgeToken DynamicGraph::add(unsigned v, unsigned u) {
        if (v == u) {
            return EdgeToken(nullptr);
        }
        unsigned n = size - 1;
        auto* edge = new Edge(n, v, u);
        if (!is_connected(v, u)) {
            edge->add_tree_edge(forests[n].link(v, u));
        }
        forests[n].increment_edges(v);
        forests[n].increment_edges(u);
        edge->subscribe(adjLists[n][v]->add(u, edge), adjLists[n][u]->add(v, edge));
        return EdgeToken(edge);
    }

    void DynamicGraph::remove(EdgeToken&& edge_token) {
        Edge* link = edge_token.edge;
        edge_token.edge = nullptr;
        if (link == nullptr) {
            return;
        }

        unsigned v = link->from();
        unsigned u = link->to();
        bool complex_deletion = link->is_tree_edge();
        unsigned level = link->level();

        if (complex_deletion) {
            for (unsigned i = 0; i <= size - level - 1; i++){
                forests[size - i - 1].cut(std::move(link->tree_edges[i]));
            }
        }

        forests[level].decrement_edges(v);
        forests[level].decrement_edges(u);

        delete link;

        if (complex_deletion) {
            for (unsigned i = level; i < size; i++){
                // find new connection
                // to do that choose the lesser component
                if(forests[i].size(v) > forests[i].size(u)){
                    std::swap(v, u);
                }
                // and iterate over good vertices until success
                // propagating all tree edges of smallest component
                Edge* replacement = nullptr;
                Iterator it = forests[i].iterator(v);
                while(it.hasNext()){
                    unsigned w = *it;
                    ListIterator lit = adjLists[i][w]->iterator();
                    while(lit.hasNext()){
                        List* l = *(lit++);
                        Edge* e = l->e();
                        unsigned up = l->vertex();
                        if (e->is_tree_edge()) {
                            downgrade(e);
                        } else {
                            if (replacement != nullptr) {
                                break;
                            }
                            if (is_connected(up, u)) {
                                replacement = e;
                            } else {
                                downgrade(e);
                            }
                        }
                    }
                    ++it;
                }

                if (replacement != nullptr) {
                    for (unsigned j = size - 1; j >= i; j--){
                        replacement->add_tree_edge(forests[j].link(replacement->v, replacement->u));
                    }
                    break;
                }
            }
        }
    }

    void DynamicGraph::downgrade(Edge* e){
        unsigned v = e->from();
        unsigned w = e->to();
        unsigned lvl = e->lvl--;
        e->removeLinks();
        e->subscribe(adjLists[lvl - 1][w]->add(v, e), adjLists[lvl - 1][v]->add(w, e));
        forests[lvl].decrement_edges(w);
        forests[lvl].decrement_edges(v);
        forests[lvl - 1].increment_edges(w);
        forests[lvl - 1].increment_edges(v);
        if (e->is_tree_edge()) {
            e->add_tree_edge(forests[lvl - 1].link(v, w));
        }
    }

    bool DynamicGraph::is_connected(unsigned v, unsigned u) {
        return forests[forests.size() - 1].is_connected(v, u);
    }

    bool DynamicGraph::is_connected() {
        return forests[forests.size() - 1].is_connected();
    }

    std::string DynamicGraph::str() {
        std::string str;
        for(unsigned i = 0; i < size; i++){
            str += "level " + std::to_string(i) + ": \n";
            str += forests[i].str() + "\n";
        }
        return str;
    }

    unsigned DynamicGraph::degree(unsigned v) {
        unsigned sum = 0;
        for (unsigned i = 0; i < size; i++) {
            sum += forests[i].degree(v);
        }
        return sum;
    }

    unsigned DynamicGraph::component_size(unsigned v) {
        return forests[forests.size() - 1].component_size(v);
    }

    List* List::add(unsigned v, Edge* edge) {
        List* new_list = new List(v, edge, prev, this);
        prev->next = new_list;
        prev = new_list;
        return new_list;
    }

    List::~List() {
        next->prev = prev;
        prev->next = next;
    }

    List::List(unsigned u, Edge* edge, List* prev, List* next) :edge(edge), u(u), prev(prev), next(next){}

    List::List() :edge(nullptr) {
        next = this;
        prev = this;
        u = std::numeric_limits<unsigned>::max();
    }

    ListIterator List::iterator() {
        return ListIterator(next);
    }

    unsigned List::vertex() {
        return u;
    }

    Edge* List::e() {
        return edge;
    }

    Edge::Edge(unsigned lvl, unsigned v, unsigned u) : lvl(lvl), v(v), u(u) {}

    void Edge::subscribe(List* first, List* second) {
        first_link = first;
        second_link = second;
    }

    unsigned Edge::level() {
        return lvl;
    }

    void Edge::removeLinks() {
        delete first_link;
        delete second_link;
        first_link = nullptr;
        second_link = nullptr;
    }

    unsigned Edge::from() {
        return v;
    }

    unsigned Edge::to() {
        return u;
    }

    void Edge::add_tree_edge(TreeEdge&& edge) {
        tree_edges.push_back(std::move(edge));
    }

    bool Edge::is_tree_edge() {
        return !tree_edges.empty();
    }

    Edge::~Edge() {
        removeLinks();
    }

    ListIterator::ListIterator(List* list) :list(list) {}

    ListIterator ListIterator::operator++(int) {
        ListIterator state(list);
        list = list->next;
        return state;
    }

    List* ListIterator::operator*() {
        return list;
    }

    bool ListIterator::hasNext() {
        return list->edge != nullptr;
    }

    EdgeToken::EdgeToken(Edge* edge) :edge(edge){}

    EdgeToken::EdgeToken(EdgeToken&& e) noexcept :edge(e.edge){
        e.edge = nullptr;
    }

    EdgeToken& EdgeToken::operator=(EdgeToken&& other) noexcept {
        edge = other.edge;
        other.edge = nullptr;
        return *this;
    }

    EdgeToken::EdgeToken() :edge(nullptr){}

    bool EdgeToken::moved() {
        return edge == nullptr;
    }
}
