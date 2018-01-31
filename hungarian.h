#ifndef HUNGARIAN_H
#define HUNGARIAN_H
#include <vector>
#include <climits>
#include <algorithm>
#include <iostream>

using std::vector;

enum class Mask{clean, starred, primed};

class Hungarian
{
public:
    Hungarian(const vector<vector<int> >& input_matrix, int rows, int cols);
    void RunHungarian();
    /* Row Reduction.
     * For each row of the matrix, find the smallest element
     * and subtract it from every element in its row.
     * Go to Step 2.
    */
    void step_one();

    /*Find a zero (Z) in the resulting matrix.
     * If there is no starred zero in its row or column, star Z.
     * Repeat for each element in the matrix.
     * Go to Step 3.
    */
    void step_two();

    /*Cover each column containing a starred zero.
     * If n columns are covered,
     * the starred zeros describe a complete set of unique assignments.
     * In this case, Go to DONE,
     * otherwise, Go to Step 4.
    */
    void step_three();

    /*Find a non-covered zero and prime it.
     * If there's no starred zero in the row of this zero, go to step 5
     * Otherwise, Cover the row and uncover the column.
     * Repeat until no uncovered zero left, go to step 6.
    */
    void step_four();

    /* This step is very similar to the augmenting path algorithm (for solving the maximal matching problem)
     * Basically we find the smallest non-covered value an subtract this from all non-covered value.
     * From step 4 we got a primed zero, find a starred zero in the column, then find a primed zero in the new row
     * Repeat until we arrive at a primed zero with no stars in the column.
     * Unstar each starred zero of the series, star each primed zero of the series,
     * erase all primes and uncover every line in the matrix.
     * Return to Step 3.
    */
    void step_five();

    /*Column Reduction
     * Add the value found in Step 4 to every element of each covered row,
     * and subtract it from every element of each uncovered column.
     * Return to Step 4 without altering any stars, primes, or covered lines.
    */
    void step_six();

    /*Done
     * Do whatever is needed after you get the result.
    */
    void step_seven();

    void find_a_nonCovered_zero(int & row, int & col);
    void find_smallest(int & min_value);
    bool star_in_row(int row);
    void find_star_in_row(int row, int & col);
    void find_star_in_col(int & row, int col);
    void find_prime_in_row(int row, int & col);
    void clearCovers();
    void erasePrimes();
    void augmentPath();

private:
    int m_rows;
    int m_cols;
    int m_step;
    int path_row_0;
    int path_col_0;
    int pathCount;

    vector<vector<int> > m_costMatrix;
    vector<vector<Mask> > m_maskMatrix;
    vector<bool> m_rowCover;
    vector<bool> m_colCover;
    vector< vector<int> > m_paths;
};

#endif // HUNGARIAN_H
