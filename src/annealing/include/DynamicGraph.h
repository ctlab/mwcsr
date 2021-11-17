#ifndef DGRAPH_DYNAMICGRAPH_H
#define DGRAPH_DYNAMICGRAPH_H

#include "EulerTourForest.h"

namespace {
    using std::vector;
}

namespace dgraph {
    class List;
    class ListIterator;
    class DynamicGraph;

    class Edge {
        unsigned lvl;
        unsigned v;
        unsigned u;
        List* first_link;
        List* second_link;
        std::vector<TreeEdge> tree_edges;
        void subscribe(List*, List*);
        void removeLinks();
        void add_tree_edge(TreeEdge&&);
    public:
        explicit Edge(unsigned, unsigned, unsigned);
        ~Edge();

        unsigned from();
        unsigned to();
        unsigned level();
        bool is_tree_edge();

        friend class DynamicGraph;
    };

    class EdgeToken {
        Edge* edge;
        explicit EdgeToken(Edge*);
    public:
        EdgeToken();
        EdgeToken(const EdgeToken&) = delete;
        EdgeToken& operator=(const EdgeToken&) = delete;
        EdgeToken& operator=(EdgeToken&&) noexcept;
        EdgeToken(EdgeToken&&) noexcept;
        ~EdgeToken() = default;

        bool moved();

        friend class DynamicGraph;
    };

    class DynamicGraph {
        unsigned n;
        unsigned size;
        vector<EulerTourForest> forests;
        vector<vector<List*>> adjLists;
        void downgrade(Edge* e);
    public:
        explicit DynamicGraph(unsigned n);
        DynamicGraph(const DynamicGraph&) = delete;
        DynamicGraph&operator=(const DynamicGraph&) = delete;
        ~DynamicGraph();

        EdgeToken add(unsigned v, unsigned u);
        void remove(EdgeToken&&);
        bool is_connected(unsigned v, unsigned u);
        bool is_connected();
        std::string str();
        unsigned degree(unsigned v);
        unsigned component_size(unsigned v);
    };

    class List {
        Edge* edge;
        unsigned u;
        List* prev;
        List* next;
        List(unsigned, Edge*, List*, List*);
    public:
        List();

        List* add(unsigned , Edge*);
        ListIterator iterator();
        unsigned vertex();
        Edge* e();
        ~List();

        friend class Edge;
        friend class ListIterator;
    };

    class ListIterator {
        List* list;
    public:
        explicit ListIterator(List*);
        ListIterator operator++(int);
        List* operator*();
        bool hasNext();
    };
}

#endif //DGRAPH_DYNAMICGRAPH_H
