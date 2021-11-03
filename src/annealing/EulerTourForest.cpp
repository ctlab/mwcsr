#include "include/EulerTourForest.h"
#include <utility>
#include <list>

namespace dgraph {

    void Entry::splay() {
        while (parent != nullptr) {
            Entry* grandpa = parent->parent;
            bool is_left = parent->left == this;
            if (grandpa != nullptr) {
                bool p_is_left = grandpa->left == parent;
                if (is_left == p_is_left) {
                    grandpa->rotate(p_is_left);
                    parent->rotate(is_left);
                } else {
                    parent->rotate(is_left);
                    grandpa->rotate(p_is_left);
                }
            } else {
                parent->rotate(is_left);
            }
        }
    }

    void Entry::remove() {
        splay();
        if (left != nullptr) {
            left->parent = nullptr;
        }
        if (right != nullptr) {
            right->parent = nullptr;
        }
        if (left == nullptr || right == nullptr){
            return;
        }
        merge(left, right);
    }

    void Entry::rotate(bool left_rotate){
        Entry* child = nullptr;
        if(left_rotate) {
            child = left;
            left = child->right;
            if (left != nullptr) {
                left->parent = this;
            }
            child->right = this;
        } else {
            child = right;
            right = child->left;
            if (right != nullptr) {
                right->parent = this;
            }
            child->left = this;
        }
        if (parent != nullptr) {
            if (this == parent->left){
                parent->left = child;
            } else {
                parent->right = child;
            }
        }
        child->parent = parent;
        parent = child;
        recalc();
        child->recalc();
        if (parent != nullptr){
            parent->recalc();
        }
    }

    Entry* merge(Entry* l, Entry* r) {
        if (l == nullptr) {
            return r;
        }
        if (r == nullptr) {
            return l;
        }
        r = find_root(r);
        l = find_root(l)->rightmost();

        l->splay();
        l->right = r;
        r->parent = l;
        l->recalc();
        return l;
    }

    Entry* find_root(Entry* e) {
        while (e->parent != nullptr) e = e->parent;
        return e;
    }

    Entry::Entry(unsigned v, Entry* l, Entry* r, Entry* p) : left(l), right(r), parent(p), v(v),
                                                             size(1), edges(0), good(false) {}

    Entry* Entry::succ() {
        Entry* curr = this;
        if(right == nullptr){
            while (curr->parent != nullptr && curr == curr->parent->right) curr = curr->parent;
            if (curr->parent == nullptr){
                return nullptr;
            }
            return curr->parent;
        }
        curr = right;
        curr = curr->leftmost();
        return curr;
    }

    std::pair<Entry*, Entry*> split(Entry* e, bool keep_in_left) {
        e->splay();
        Entry* left;
        Entry* right;
        if (keep_in_left) {
            left = e;
            right = e->right;
            e->right = nullptr;
            left->recalc();
            if (right != nullptr) {
                right->recalc();
                right->parent = nullptr;
            }
        } else {
            left = e->left;
            right = e;
            e->left = nullptr;
            right->recalc();
            if (left != nullptr) {
                left->recalc();
                left->parent = nullptr;
            }
        }
        return std::make_pair(left, right);
    }

    void Entry::recalc() {
        size = 1;
        good = edges > 0;
        if(right != nullptr){
            size += right->size;
            good |= right->good;
        }
        if(left != nullptr){
            size += left->size;
            good |= left->good;
        }
    }

    EulerTourForest::EulerTourForest(unsigned n) : n(n), any_root(nullptr) {
        for (unsigned i = 0; i < n; i++) {
            auto* vertex = new Entry(i);
            any.push_back(vertex);
        }
    }

    EulerTourForest::EulerTourForest(EulerTourForest&& forest) noexcept :n(forest.n), any(std::move(forest.any)),
                                                                         any_root(forest.any_root) {
        forest.n = 0;
    }

    EulerTourForest::~EulerTourForest() {
        std::vector<bool> vis(n, false);
        std::list<Entry*> entries;
        for (int i = 0; i < n; i++) {
            if (vis[i]) {
                continue;
            }
            vis[i] = true;
            Entry* e = find_root(any[i])->leftmost();
            while (e != nullptr){
                vis[e->v] = true;
                entries.push_back(e);
                e = e->succ();
            }
        }
        for (Entry* e : entries) {
            delete e;
        }
    }

    Entry* EulerTourForest::make_root(unsigned v) {
        Entry* e = any[v];
        auto cut = split(e, false);
        return merge(cut.second, cut.first);
    }

    Entry* EulerTourForest::expand(unsigned v) {
        Entry* e = make_root(v);
        if (e->size == 1){
            return e;
        }
        auto new_node = new Entry(v);
        merge(e, new_node);
        return new_node;
    }

    TreeEdge EulerTourForest::link(unsigned v, unsigned u) {
        Entry* l = expand(v);
        Entry* r = expand(u);
        any_root = merge(l, r);
        return {l, r};
    }

    void EulerTourForest::cut(Entry* first, Entry* last) {
        any_root = nullptr;
        auto first_cut = split(first, true);
        bool right_ordered = first_cut.second != nullptr && find_root(first_cut.second) == find_root(last);
        auto second_cut = split(last, true);
        if (!right_ordered) {
            std::swap(first_cut, second_cut);
        }
        Entry* to_remove = first_cut.first->rightmost();
        if (to_remove->is_singleton()) {
            if (second_cut.second != nullptr) {
                change_any(second_cut.second->leftmost());
                delete to_remove;
            }
        } else {
            merge(to_remove, second_cut.second);
            Entry* next = to_remove->succ();
            if (next == nullptr) {
                cutoff(to_remove);
            } else {
                cutoff(to_remove, next);
            }
        }
        cutoff(second_cut.first->rightmost());
    }

    void EulerTourForest::cutoff(Entry* e, Entry* replacement) {
        if (e->is_singleton()) {
            return;
        }
        if (any[e->v] == e){
            if (replacement == nullptr) {
                change_any(find_root(e)->leftmost());
            } else {
                change_any(replacement);
            }
        }
        e->remove();
        delete e;
    }

    void EulerTourForest::change_any(Entry* e) {
        unsigned edges = any[e->v]->edges;
        unsigned v = e->v;
        change_edges(v, 0);
        any[v] = e;
        change_edges(v, edges);
    }

    bool EulerTourForest::is_connected() {
        return any_root != nullptr && (int)any_root->size == 2 * (n - 1);
    }

    bool EulerTourForest::is_connected(unsigned v, unsigned u) {
        if (is_connected()) {
            return true;
        }
        return find_root(any[v]) == find_root(any[u]);
    }

    void EulerTourForest::increment_edges(unsigned v) {
        Entry* curr = any[v];
        ++curr->edges;
        if (curr->edges == 1) {
            curr->good = true;
            repair_edges_number(curr->parent);
        }
    }

    void EulerTourForest::decrement_edges(unsigned v) {
        Entry* curr = any[v];
        --curr->edges;
        if (curr->edges == 0) {
            repair_edges_number(curr);
        }
    }

    void EulerTourForest::change_edges(unsigned v, unsigned n) {
        Entry* curr = any[v];
        curr->edges = n;
        repair_edges_number(curr);
    }

    void EulerTourForest::repair_edges_number(Entry* curr){
        while (curr != nullptr) {
            bool good = curr->edges > 0;
            if (curr->left != nullptr) {
                good |= curr->left->good;
            }
            if (curr->right != nullptr) {
                good |= curr->right->good;
            }
            if (good != curr->good) {
                curr->good = good;
                curr = curr->parent;
            } else {
                return;
            }
        }
    }

    unsigned EulerTourForest::size(unsigned v) {
        return find_root(any[v])->size;
    }

    Iterator EulerTourForest::iterator(unsigned v){
        return any[v]->iterator();
    }

    std::string EulerTourForest::str() {
        std::string str;
        std::vector<bool> vis(n, false);
        for (int i = 0; i < n; i++) {
            Entry* curr = find_root(any[i]);
            if(!vis[curr->vertex()]){
                vis[curr->vertex()] = true;
                str += curr->str() + "\n";
            }
        }
        str += "edges: \n";
        for (int i = 0; i < n; i++) {
            str += std::to_string(any[i]->edges) + " ";
        }
        str += "\n";
        return str;
    }

    void EulerTourForest::cut(TreeEdge&& edge) {
        if (edge.edge != nullptr) {
            cut(edge.edge, edge.twin);
        }
    }

    unsigned EulerTourForest::degree(unsigned v) {
        return any[v]->edges;
    }

    unsigned EulerTourForest::component_size(unsigned v) {
        unsigned nodes = size(v);
        if (nodes == 1) {
            return 1;
        }
        return nodes / 2 + 1;
    }

    Iterator::Iterator(Entry* entry) :entry(entry){}

    Iterator& Iterator::operator++() {
        if (entry->right != nullptr && entry->right->good){
            entry = entry->right;
            while (true){
                if (entry->left != nullptr && entry->left->good){
                    entry = entry->left;
                    continue;
                }
                if(entry->edges > 0){
                    return *this;
                }
                entry = entry->right;
            }
        }
        while (true) {
            if (entry->parent == nullptr) {
                entry = nullptr;
                return *this;
            }
            if (entry->parent->right != nullptr && entry->parent->right == entry){
                entry = entry->parent;
                continue;
            } else {
                entry = entry->parent;
                break;
            }
        }
        if (entry->edges > 0){
            return *this;
        }
        return ++(*this);
    }

    unsigned Iterator::operator*() {
        return entry->v;
    }

    bool Iterator::hasNext() {
        return entry != nullptr;
    }

    Iterator Entry::iterator() {
        Entry* curr = find_root(this)->leftmost();
        Iterator iterator(curr);
        if(!curr->good) {
            ++iterator;
        }
        return iterator;
    }

    unsigned Entry::vertex() {
        return v;
    }

    std::string Entry::str() {
        std::string str;
        Entry* e = leftmost();
        while(e != nullptr){
            str += std::to_string(e->v);
            e = e->succ();
        }
        return str;
    }

    Entry* Entry::leftmost() {
        Entry* curr = this;
        while (curr->left != nullptr) curr = curr->left;
        return curr;
    }

    Entry* Entry::rightmost() {
        Entry* curr = this;
        while (curr->right != nullptr) curr = curr->right;
        return curr;
    }

    bool Entry::is_singleton() {
        return parent == nullptr && left == nullptr && right == nullptr;
    }

    TreeEdge::TreeEdge(Entry* e, Entry* t) :edge(e), twin(t) {}

    TreeEdge::TreeEdge(TreeEdge&& edge) noexcept :edge(edge.edge), twin(edge.twin){
        edge.twin = nullptr;
        edge.edge = nullptr;
    }
}
