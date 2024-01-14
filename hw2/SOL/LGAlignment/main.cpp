#include <string>
#include <iostream>
#include "argh.h"
#include <sstream>
#include <fstream>
#include <climits>

using namespace std;

int main(int argc, char **argv) {

    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    bool global_alignment;

    if(cmdl["-g"])
        global_alignment = true;
    else if( cmdl["-l"])
        global_alignment = false;
    else {
        cout << "No valid flag for alignment type. Terminating..." << endl;
        exit(1);
    }


    string pattern_file, text_file, output_file;
    cmdl("-p") >> pattern_file;
    cmdl("-t") >> text_file;
    cmdl("-o") >> output_file;

    string scoring;
    cmdl("-s") >> scoring;

    string scores[3];
    int index = 0;

    string word = "";
    for (auto x : scoring)
    {
        if (x == '|')
        {
            scores[index++] = word;
            word = "";
        }
        else {
            word = word + x;
        }
    }

    scores[index] = word;


    int match = stoi(scores[0]), mismatch_penalty = stoi(scores[1]), gap_penalty = stoi(scores[2]);

    ifstream input_pattern(pattern_file);
    ifstream input_text(text_file);

    string seq1; // pattern
    getline(input_pattern, seq1);
    input_pattern.close();

    string seq2; // text
    getline(input_text, seq2);
    input_text.close();

    const int dim1 = seq1.length() + 1;
    const int dim2 = seq2.length() + 1;

    int alignment_matrix[dim1][dim2];
    int match_matrix[dim1 - 1][dim2 - 1];

    int score;

    if(global_alignment){
        // Setting alignment matrix all zeros!
        for(int i = 0; i < dim1; i++){

            for(int j = 0; j < dim2; j++){
                alignment_matrix[i][j] = 0;
            }
        }

        // Setting first row & first col with gap penalties
        for(int i = 0; i < dim1; i++){
            alignment_matrix[i][0] = i * gap_penalty;
        }

        for(int i = 0; i < dim2; i++){
            alignment_matrix[0][i] = i * gap_penalty;
        }

        for(int i = 0; i < dim1 - 1; i++){

            for(int j = 0; j < dim2 - 1; j++){
                match_matrix[i][j] = (seq1.at(i) == seq2.at(j)) ? match:mismatch_penalty;
            }
        }

        // Starting from index 1, filling the matrix with the max score of i-1, j-1 + match points or gaps from i-1, j or i, j-1
        for(int i = 1; i < dim1; i++){

            for(int j = 1; j < dim2; j++){

                int max_score = alignment_matrix[i - 1][j] + gap_penalty;
                max_score = max(max_score, alignment_matrix[i][j - 1] + gap_penalty);
                max_score = max(max_score, alignment_matrix[i - 1][j - 1] + match_matrix[i - 1][j - 1]);

                alignment_matrix[i][j] = max_score;
            }
        }

        score = alignment_matrix[dim1 - 1][dim2 - 1];

        ofstream output(output_file);

        for(size_t i = 0; i < seq2.length(); i++ ){
            if( i == 0) {
                output.width(10);
                output <<  " ";
            }
            output.width(7);
            output << seq2.at(i) << " \n"[ i == seq2.length() - 1];
        }

        for(int i = 0; i < dim1; i++){
            if( i == 0){
                output <<  "  ";
            }
            if( i > 0) {
                output << seq1.at(i - 1) << " ";
            }

            for(int j = 0; j < dim2; j++){

                output.width(7);
                output << alignment_matrix[i][j] << " \n"[ j == dim2 - 1];
            }
        }

        output.close();

        // Alignment
        string row_aligned = "", col_aligned = "";

        int row = dim1 - 1, col = dim2 - 1;

        while( row > 0 || col > 0){
            if( row == 0){
                row_aligned = "-" + row_aligned;
                col_aligned = seq2[col - 1] + col_aligned;
                col--;
            }
            else if(col == 0) {
                row_aligned = seq1[row - 1] + row_aligned;
                col_aligned = "-" + col_aligned;
                row--;
            }
            else {
                if( alignment_matrix[row][col] == alignment_matrix[row - 1][col - 1] + match_matrix[row - 1][col - 1]){
                    row_aligned = seq1[row - 1] + row_aligned;
                    col_aligned = seq2[col - 1] + col_aligned;
                    row--;
                    col--;
                }
                else if( alignment_matrix[row][col] == alignment_matrix[row][col - 1] + gap_penalty )
                {
                    row_aligned = "-" + row_aligned;
                    col_aligned = seq2[col - 1] + col_aligned;
                    col--;
                }
                else {

                    row_aligned = seq1[row - 1] + row_aligned;
                    col_aligned = "-" + col_aligned;
                    row--;
                }
            }
        }

        cout << row_aligned << endl;
        cout << col_aligned << endl;

    }
    else {

        for(int i = 0; i < dim1 - 1; i++){

            for(int j = 0; j < dim2 - 1; j++){
                match_matrix[i][j] = (seq1.at(i) == seq2.at(j)) ? match:mismatch_penalty;
            }
        }

        int ind_score;
        for(int i = 0; i < dim1; i++){

            ind_score = (i * gap_penalty < 0) ? 0:i * gap_penalty;
            alignment_matrix[i][0] = ind_score;
        }

        for(int i = 0; i < dim2; i++){
            ind_score = (i * gap_penalty < 0) ? 0:i * gap_penalty;
            alignment_matrix[0][i] = ind_score;
        }

        for(int i = 1; i < dim1; i++){

            for(int j = 1; j < dim2; j++){

                int max_score = alignment_matrix[i - 1][j] + gap_penalty;
                max_score = max(max_score, alignment_matrix[i][j - 1] + gap_penalty);
                max_score = max(max_score, alignment_matrix[i - 1][j - 1] + match_matrix[i - 1][j - 1]);

                ind_score = (max_score < 0) ? 0 : max_score;

                alignment_matrix[i][j] = ind_score;
            }
        }

        ofstream output(output_file);

        for(size_t i = 0; i < seq2.length(); i++ ){
            if( i == 0) {
                output.width(10);
                output <<  " ";
            }
            output.width(7);
            output << seq2.at(i) << " \n"[ i == seq2.length() - 1];
        }

        for(int i = 0; i < dim1; i++){
            if( i == 0){
                output <<  "  ";
            }
            if( i > 0) {
                output << seq1.at(i - 1) << " ";
            }

            for(int j = 0; j < dim2; j++){

                output.width(7);
                output << alignment_matrix[i][j] << " \n"[ j == dim2 - 1];
            }
        }

        output.close();

        // Alignment

        // Finding maximum match score in local alignment matrix
        int max_score = INT_MIN;
        int max_row, max_col;

        for(int i = 0; i < dim1; i++){
            for(int j = 0; j < dim2; j++){
                if(alignment_matrix[i][j] > max_score){
                    max_score = alignment_matrix[i][j];
                    max_row = i;
                    max_col = j;
                }
            }
        }


        score = max_score;

        string row_aligned = "", col_aligned = "";
        int row = max_row, col = max_col;
        /*while(row > 0 && col > 0 && alignment_matrix[row][col] != 0){
            if( alignment_matrix[row][col] == alignment_matrix[row - 1][col - 1] + match_matrix[row - 1][col - 1]){
                    row_aligned = seq1[row - 1] + row_aligned;
                    col_aligned = seq2[col - 1] + col_aligned;
                    row--;
                    col--;
            }
        }*/

        while( row > 0 && col > 0 && alignment_matrix[row][col] != 0){

                if( alignment_matrix[row][col] == alignment_matrix[row - 1][col - 1] + match_matrix[row - 1][col - 1]){
                    row_aligned = seq1[row - 1] + row_aligned;
                    col_aligned = seq2[col - 1] + col_aligned;
                    row--;
                    col--;
                }
                else if( alignment_matrix[row][col] == alignment_matrix[row][col - 1] + gap_penalty )
                {
                    row_aligned = "-" + row_aligned;
                    col_aligned = seq2[col - 1] + col_aligned;
                    col--;
                }
                else {

                    row_aligned = seq1[row - 1] + row_aligned;
                    col_aligned = "-" + col_aligned;
                    row--;
                }
        }


        cout << row_aligned << endl;
        cout << col_aligned << endl;
    }

    cout << "Score=" << score << endl;
    return 0;
}
