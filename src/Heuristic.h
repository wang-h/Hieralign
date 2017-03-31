#ifndef HEURISTIC_H
#define HEURISTIC_H

#include "Matrix.h"
const vector<int> offset{-1, 1};
inline bool AddVerticalNeighbors(const unsigned &i, const unsigned &j, const unsigned &m, Matrix *matrix,
                                 set<unsigned> *F_aligned, Alignment *sure)
{
    bool added = false;
    for (auto &p : offset)
    {
        const unsigned x = i + p;
        if ((1 < i) && (i < m - 1) && F_aligned->find(x) == F_aligned->end() && (*matrix)(x, j) == 1)
        {
            sure->push_back(make_pair(j, x));
            (*matrix)(x, j) = 2;
            F_aligned->insert(x);
            added = true;
        }
    }
    return added;
}
inline bool AddHorizontalNeighbors(const unsigned &i, const unsigned &j, const unsigned &n, Matrix *matrix,
                                   set<unsigned> *E_aligned, Alignment *sure)
{
    bool added = false;
    for (auto &p : offset)
    {
        const unsigned y = j + p;
        if ((1 < j) && (j < n - 1) && E_aligned->find(y) == E_aligned->end() && (*matrix)(i, y) == 1)
        {
            sure->push_back(make_pair(y, i));
            (*matrix)(i, y) = 2;
            E_aligned->insert(y);
            added = true;
        }
    }
    return added;
}
inline bool AddDiagonalNeighbors(const unsigned &i, const unsigned &j, const unsigned &m, const unsigned &n, Matrix *matrix,
                                 set<unsigned> *F_aligned, set<unsigned> *E_aligned, Alignment *sure)
{
    bool added = false;
    for (auto &p : offset)
    {
        for (auto &q : offset)
        {
            const unsigned x = i + p;
            const unsigned y = j + q;
            if ((1 < j) && (j < n - 1) && (1 < i) && (i < m - 1))
            {
                if (F_aligned->find(i) == F_aligned->end() || E_aligned->find(j) == E_aligned->end())
                {
                    if ((*matrix)(x, y) == 1)
                    {
                        sure->push_back(make_pair(y, x));
                        (*matrix)(x, y) = 2;
                        F_aligned->insert(x);
                        E_aligned->insert(y);
                        added = true;
                    }
                }
            }
        }
    }
    return added;
}

inline void AddFinalAndNeighbors(const unsigned &m, const unsigned &n, Matrix *matrix,
                                 set<unsigned> *F_aligned, set<unsigned> *E_aligned, Alignment *sure)
{
    for (unsigned i = 0; i < m; i++)
    {
        for (unsigned j = 0; j < n; j++)
        {
            if (F_aligned->find(i) == F_aligned->end() && E_aligned->find(j) == E_aligned->end())
            {
                if ((*matrix)(i, j) == 1)
                {
                    sure->push_back(make_pair(j, i));
                    F_aligned->insert(i);
                    E_aligned->insert(j);
                }
            }
        }
    }
}
void Heuristic(const unsigned &m, const unsigned &n, set<unsigned> *F_aligned, set<unsigned> *E_aligned,
               Alignment *sure, Matrix *matrix)
{
    bool added = false;
    while (true)
    {
        added = false;
        for (unsigned i = 0; i < m; i++)
        {
            for (unsigned j = 0; j < n; j++)
            {
                if ((*matrix)(i, j) == 2)
                {
                    if (AddVerticalNeighbors(i, j, m, matrix, F_aligned, sure) ||
                        AddHorizontalNeighbors(i, j, n, matrix, F_aligned, sure) ||
                        AddDiagonalNeighbors(i, j, m, n, matrix, F_aligned, E_aligned, sure))
                        added = true;
                }
            }
        }
        if (!added)
            break;
    }
    AddFinalAndNeighbors(m, n, matrix, F_aligned, E_aligned, sure);
}

void SymmetrizeAlignment(const AlignmentList &forward_list, const AlignmentList &backward_list,
                         AlignmentList *symmetric_list)
{
    CHECK(forward_list.size() == backward_list.size(), "#ERROR, two alignment lists with different size!");
    symmetric_list->resize(forward_list.size());
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < symmetric_list->size(); i++)
    {

        const Alignment &forward_links = forward_list[i];
        const Alignment &backward_links = backward_list[i];
        const unsigned n = static_cast<unsigned>(forward_links.size());
        const unsigned m = static_cast<unsigned>(backward_links.size());
        Matrix align_matrix(m, n);
        Alignment sure;
        set<unsigned> F_aligned, E_aligned;
        for (auto &it : forward_links)
        {
            const unsigned int y = it.first;
            const unsigned int x = it.second;
            if (x < m)
                align_matrix(x, y) = 1;
        }
        for (auto &it : backward_links)
        {
            const unsigned int x = it.first;
            const unsigned int y = it.second;
            if (y < n)
            {
                if (align_matrix(x, y) > 0)
                {
                    sure.push_back(make_pair(y, x));
                    F_aligned.insert(x);
                    E_aligned.insert(y);
                }
                align_matrix(x, y) += 1;
            }
        }
        //align_matrix.PrintMatrix();
        Heuristic(m, n, &F_aligned, &E_aligned, &sure, &align_matrix);

        (*symmetric_list)[i] = sure;
    }
}

#endif // HEURISTIC_H
