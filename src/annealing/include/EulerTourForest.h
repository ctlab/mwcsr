#ifndef DGRAPH_EULERTOURTREE_H
#define DGRAPH_EULERTOURTREE_H

#include <vector>
#include <string>

namespace dgraph {
    class Iterator;
    class EulerTourForest;

    class Entry {
        Entry* left;
        Entry* right;
        Entry* parent;
        unsigned v;
        unsigned size;
        unsigned edges;
        bool good;

        explicit Entry(unsigned , Entry* = nullptr, Entry* = nullptr, Entry* = nullptr);

        void splay();
        void rotate(bool);
        void remove();
        Entry* succ();
        Entry* leftmost();
        Entry* rightmost();
        void recalc();
        Iterator iterator();
        bool is_singleton();
        std::string str();
        unsigned depth(unsigned);

        friend Entry* merge(Entry*, Entry*);
        friend std::pair<Entry*, Entry*> split(Entry*, bool);
        friend Entry* find_root(Entry* e);

        friend class EulerTourForest;
        friend class Iterator;

    public:
        unsigned vertex();
    };

    class Iterator {
        Entry* entry;
    public:
        explicit Iterator(Entry*);
        Iterator& operator++();
        unsigned operator*();
        bool hasNext();
    };

    class TreeEdge {
        Entry* edge;
        Entry* twin;
        TreeEdge(Entry*, Entry*);
    public:
        TreeEdge(TreeEdge&&) noexcept;
        ~TreeEdge() = default;

        friend class EulerTourForest;
    };

    class EulerTourForest {
        int n;
        std::vector<Entry*> any;
        Entry* any_root;
        Entry* make_root(unsigned v);
        Entry* expand(unsigned v);
        void change_any(Entry* e);
        void cutoff(Entry* e, Entry* replacement = nullptr);
        void cut(Entry*, Entry*);
        void repair_edges_number(Entry*);

    public:
        explicit EulerTourForest(unsigned);
        EulerTourForest(const EulerTourForest&) = delete;
        EulerTourForest& operator=(const EulerTourForest&) = delete;
        EulerTourForest(EulerTourForest&&) noexcept;
        ~EulerTourForest();

        bool is_connected(unsigned v, unsigned u);
        bool is_connected();
        TreeEdge link(unsigned v, unsigned u);
        void cut(TreeEdge&&);
        void increment_edges(unsigned v);
        void decrement_edges(unsigned v);
        void change_edges(unsigned v, unsigned n);
        unsigned size(unsigned v);
        Iterator iterator(unsigned v);
        std::string str();
        unsigned degree(unsigned v);
        unsigned component_size(unsigned v);
    };
}

#endif //DGRAPH_EULERTOURTREE_H
