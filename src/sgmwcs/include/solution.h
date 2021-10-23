#ifndef MWCSR_SOLUTION_H
#define MWCSR_SOLUTION_H

#include <vector>
#include <cstddef>

namespace relax {

class Solution {
    std::vector<size_t> edges;
    double obj;
public:
    Solution();
    Solution(std::vector<size_t> elements, double obj);
    Solution(const Solution&) = default;
    Solution& operator=(const Solution&) = default;
    Solution(Solution&&) = default;

    void add_edge(size_t e);
    void set_obj(double obj);
    void add_obj(double value);

    std::vector<size_t> solution() const;
    double objective() const;

};

}

#endif //MWCSR_SOLUTION_H
