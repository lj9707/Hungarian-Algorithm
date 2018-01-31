#include "hungarian.h"

Hungarian::Hungarian(const vector<vector<int> >& input_matrix, int rows, int cols):
    m_costMatrix(input_matrix),
    m_rows(rows),
    m_cols(cols),
    m_step(1),
    m_rowCover(vector<bool>(rows, false)),
    m_colCover(vector<bool>(cols, false)),
    m_maskMatrix(vector< vector<Mask> >(
                                       rows,
                                       vector<Mask>(cols, Mask::clean)
                                      )),
    m_paths(vector< vector<int> >(
                                    rows*rows,
                                    vector<int>(2, -2)
                                ))
{

}

void Hungarian::RunHungarian()
{
    bool done = false;
    while(!done)
    {
        switch (m_step)
        {
            case 1:
                step_one();
                break;
            case 2:
                step_two();
                break;
            case 3:
                step_three();
                break;
            case 4:
                step_four();
                break;
            case 5:
                step_five();
                break;
            case 6:
                step_six();
                break;
            case 7:
                step_seven();
                done = true;
                break;
            default:

                break;
        }
    }
}

void Hungarian::step_one()
{
    int minOfRow = INT_MAX;
    for(int r = 0; r < m_rows; ++r)
    {
        //find min of current row
        minOfRow = m_costMatrix[r][0];
        for(int c = 0; c < m_cols; ++c)
        {
            minOfRow = std::min(minOfRow, m_costMatrix[r][c]);
        }
        //subtract every element by the min
        for(int c = 0; c < m_cols; ++c)
        {
            m_costMatrix[r][c] -= minOfRow;
        }
    }
    m_step = 2;
}

void Hungarian::step_two()
{
    for(int r = 0; r < m_rows; ++r)
    {
        for(int c = 0; c < m_cols; ++c)
        {
            if(m_costMatrix[r][c] == 0
                    && !m_rowCover[r]
                    && !m_colCover[c])
            {
                m_maskMatrix[r][c] = Mask::starred;
                m_rowCover[r] = true;
                m_colCover[c] = true;
            }
        }
    }

    clearCovers();

    m_step = 3;
}

void Hungarian::step_three()
{
    int colCount = 0;
    for(int r = 0; r < m_rows; ++r)
        for(int c = 0; c < m_cols; ++c)
            if(m_maskMatrix[r][c] == Mask::starred)
                m_colCover[c] = true;
    for(int c = 0; c < m_cols; ++c)
        if(m_colCover[c])
            ++colCount;
    if(colCount >= m_cols)
        m_step = 7;
    else
        m_step = 4;
}

void Hungarian::step_four()
{
    int row = -1;
    int col = -1;

    bool done = false;
    while(!done)
    {
        find_a_nonCovered_zero(row, col);
        if(row == -1)
        {
            done = true;
            m_step = 6;
        }else
        {
            //if(m_maskMatrix[row][col] != Mask::starred)
                m_maskMatrix[row][col] = Mask::primed;
            if(star_in_row(row))
            {
                find_star_in_row(row,col);
                m_rowCover[row] = true;
                m_colCover[col] = false;
            }
            else
            {
                done = true;
                m_step = 5;
                path_row_0 = row;
                path_col_0 = col;
            }

        }
    }


}

void Hungarian::step_five()
{
    bool done = false;
    int r = -1;
    int c = -1;
    pathCount = 1;
    m_paths[pathCount - 1][0] = path_row_0;
    m_paths[pathCount - 1][1] = path_col_0;

    while(!done)
    {
        find_star_in_col(r, m_paths[pathCount-1][1] );
        if(r != -1)
        {
            ++pathCount;
            m_paths[pathCount-1][0] = r;
            m_paths[pathCount-1][1] = m_paths[pathCount-2][1];
        }
        else
        {
            done = true;
        }
        if(!done)
        {
            find_prime_in_row(m_paths[pathCount-1][0], c);
            ++pathCount;
            m_paths[pathCount-1][0] = m_paths[pathCount-2][0];
            m_paths[pathCount-1][1] = c;
        }
    }
    augmentPath();
    clearCovers();
    erasePrimes();
    m_step = 3;
}

void Hungarian::step_six()
{
    int min_value = INT_MAX;
    find_smallest(min_value);

    for(int r = 0; r < m_rows; ++r)
    {
        for(int c = 0; c < m_cols; ++c)
        {
            if(m_rowCover[r])
                m_costMatrix[r][c] += min_value;
            if(!m_colCover[c])
                m_costMatrix[r][c] -= min_value;
        }
    }
    m_step = 4;
}

void Hungarian::step_seven()
{
    for(int r = 0; r < m_rows; ++r)
    {
        for(int c = 0; c < m_cols; ++c)
        {
            if(m_maskMatrix[r][c] == Mask::starred)
            {
                std::cout << "*|";
            }
            else
            {
                std::cout << " |";
            }
            std::cout << " ";
        }
        std::cout << "\n";
    }
}

void Hungarian::find_a_nonCovered_zero(int &row, int &col)
{
    row = -1;
    col = -1;
    bool done = false;
    for(int r = 0; r < m_rows; ++r)
    {
        for(int c = 0; c < m_cols; ++c)
        {
            if(m_costMatrix[r][c] == 0
                    && (!m_rowCover[r])
                    && (!m_colCover[c]))
            {
                row = r;
                col = c;
                done = true;
            }
            if(done)
                break;
        }
        if(done)
            break;
    }
}

void Hungarian::find_smallest(int &min_value)
{
    for(int r = 0; r < m_rows; ++r)
        for(int c = 0; c < m_cols; ++c)
            if((!m_rowCover[r]) && (!m_colCover[c]))
                min_value = std::min(min_value, m_costMatrix[r][c]);
}

bool Hungarian::star_in_row(int row)
{
    for(int c = 0; c < m_cols; ++c)
        if(m_maskMatrix[row][c] == Mask::starred)
            return true;
    return false;
}

void Hungarian::find_star_in_row(int row, int &col)
{
    col = -1;
    for(int c = 0; c < m_cols; ++c)
        if(m_maskMatrix[row][c] == Mask::starred)
            col = c;
}

void Hungarian::find_star_in_col(int &row, int col)
{
    row = -1;
    for(int r = 0; r < m_rows; ++r)
        if(m_maskMatrix[r][col] == Mask::starred)
            row = r;
}

void Hungarian::find_prime_in_row(int row, int &col)
{
    for(int c = 0; c < m_cols; ++c)
        if(m_maskMatrix[row][c] == Mask::primed)
            col = c;
}

void Hungarian::clearCovers()
{
    std::fill(m_rowCover.begin(), m_rowCover.end(), false);
    std::fill(m_colCover.begin(), m_colCover.end(), false);
}

void Hungarian::erasePrimes()
{
    for(int r = 0; r < m_rows; ++r)
        for(int c = 0; c < m_cols; ++c)
            if(m_maskMatrix[r][c] == Mask::primed)
                m_maskMatrix[r][c] = Mask::clean;
}

void Hungarian::augmentPath()
{
    for(int p = 0; p < pathCount; ++p)
    {
        int p1 = m_paths[p][0];
        int p2 = m_paths[p][1];
        if(m_maskMatrix[p1][p2] == Mask::starred)
            m_maskMatrix[p1][p2] = Mask::clean;
        else
            m_maskMatrix[p1][p2] = Mask::starred;
    }
}

