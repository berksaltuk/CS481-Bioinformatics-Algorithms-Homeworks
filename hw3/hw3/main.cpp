#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <climits>
#include <tuple>

#include "argh.h"
using namespace std;

tuple<string, string> needleman_wunsch(string seq1, string seq2, int match, int mismatch_penalty, int gap_penalty) {

    const int dim1 = seq1.length() + 1;
    const int dim2 = seq2.length() + 1;

    int alignment_matrix[dim1][dim2];
    int match_matrix[dim1 - 1][dim2 - 1];


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

    return  make_tuple(row_aligned, col_aligned);
}

int calculate_edit_distance(string seq1, string seq2, int match, int mismatch_penalty, int gap_penalty){

    string row_aligned, col_aligned;

    tie(row_aligned, col_aligned) = needleman_wunsch(seq1, seq2, match, mismatch_penalty, gap_penalty);

    int distance = 0;

    for(size_t i = 0; i < row_aligned.length(); i++){
        if(row_aligned.at(i) != col_aligned.at(i)){
            distance++;
        }
    }

    return distance;
}

string merge_centers(string prev_center, string cur_center){

    string merged = "";
    size_t i = 0, j = 0;

    while(i < prev_center.length() && j < cur_center.length()){
        if(prev_center.at(i) == cur_center.at(j)){
            merged += prev_center.at(i);
            i++, j++;
        }
        else if(prev_center.at(i) == '-'){
            merged += "-";
            i++;
        }
        else if(cur_center.at(j) == '-'){
            merged += "-";
            j++;
        }
    }

    while(i < prev_center.length()){
        merged += prev_center.at(i);
        i++;
    }

    while(j < cur_center.length()){
        merged += cur_center.at(j);
        j++;
    }

    return merged;
}

string modify_pattern(string pattern, string ultimate_center, string center_i){

    size_t c_i = 0;

    for(size_t i = 0; i < ultimate_center.length(); i++){
        if(c_i < center_i.length() && ultimate_center.at(i) == '-' && pattern.at(i) != '-' && center_i.at(c_i) != '-')
        {
            pattern.insert(i, "-");
        } else {
            c_i++;
        }
    }

    while(pattern.length() < ultimate_center.length()){
        pattern += "-";
    }

    return pattern;
}


int main(int argc, char **argv){

    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    string input_fasta, output_phy;
    cmdl("-i") >> input_fasta;
    cmdl("-o") >> output_phy;

    string scoring;
    cmdl("-s") >> scoring;

    string scores[3];
    int index = 0;

    string word = "";
    for (auto x : scoring)
    {
        if (x == ':')
        {
            scores[index++] = word;
            word = "";
        }
        else {
            word = word + x;
        }
    }

    scores[index] = word;

    int match = stoi(scores[0]);
    int mismatch_penalty = stoi(scores[1]);
    int gap_penalty = stoi(scores[2]);


    int seq_count = 0;

    string linello;

    ifstream ifstream1(input_fasta);
    while( getline ( ifstream1, linello ) )
    {
        if(linello.at(0) == '>')
            seq_count++;
    }
    ifstream1.close();

    string seq_arr[seq_count];
    string seq_names[seq_count];

    for(int i = 0; i < seq_count; i++){
        seq_arr[i] = "";
    }

    ifstream ifstream2(input_fasta);
    int cur = 0, cur2 = -1;
    while( getline ( ifstream2, linello ) )
    {
        if(linello.at(0) == '>'){
            seq_names[cur++] = linello.substr(1);
            cur2++;
        }
        else{
            seq_arr[cur2] += linello;
        }
    }
    ifstream2.close();

    cout << "--- IN CASE WRITING TO OUTPUT FILE GOES WRONG: CONSOLE OUTPUTS! ---" << endl;

    cout << endl;
    int distance_matrix[seq_count][seq_count];

    for(int i = 0; i < seq_count; i++){
        for( int j = 0; j < seq_count; j++){
            distance_matrix[i][j] = calculate_edit_distance(seq_arr[i], seq_arr[j], match, mismatch_penalty, gap_penalty);
        }
    }

    cout << "DISTANCE MATRIX FOR CHOOSING THE CENTER SEQUENCE: " << endl;
    cout << endl;
    for (int i = 0; i < seq_count; i++)
       for (int j = 0; j < seq_count; j++)
          cout << distance_matrix[i][j] << " \n"[j == seq_count-1];

    cout << endl;

    int min_sum = INT_MAX;
    int row_sum = 0;
    int min_index = 0;

    for (int i = 0; i < seq_count; i++){
        for (int j = 0; j < seq_count; j++){
            row_sum += distance_matrix[i][j];
        }
        if(row_sum < min_sum){
            min_index = i;
            min_sum = row_sum;
        }
        row_sum = 0;
    }

    string center = seq_arr[min_index];

    cout << "--------------------------------------------------------" << endl;
    cout << "CENTER SEQUENCE WITH MINIMUM DISTANCE TO OTHER SEQUENCES: " << endl;
    cout << endl;

    cout << "Sequence #: " << min_index << ", Total Distance: " << min_sum << endl;
    cout << endl;

    cout << "Center Sequence: " << endl;
    cout << endl;

    cout << center << endl;
    cout << endl;


    int pair_size = seq_count - 1;
    string centers[pair_size];
    string aligned_pairs[pair_size];



    int pair_i = 0;

    for(int i = 0; i < seq_count; i++){
        if( i != min_index) {
            tie(centers[pair_i], aligned_pairs[pair_i]) = needleman_wunsch(center, seq_arr[i], match, mismatch_penalty, gap_penalty);
            pair_i++;
        }
    }

    string ultimate_center = centers[0];

    for(int i = 0; i < pair_size; i++)
    {
        if(i != 0)
            ultimate_center = merge_centers(ultimate_center, centers[i]);

        for(int k = 0; k <= i; k++){
            aligned_pairs[k] = modify_pattern(aligned_pairs[k], ultimate_center, centers[k]);
            centers[k] = ultimate_center;
        }
    }

    cout << "-------------------------------------------------------" << endl;
    cout << endl;

    cout << "ALIGNED SEQUENCES (THEY DO NOT LOOK SO GOOD SO I WISH YOU CAN SEE THE OUTPUT FILE): " << endl;
    cout << endl;

    cout << ultimate_center << endl;
    cout << endl;

    for(int i = 0; i < pair_size; i++){
        cout << aligned_pairs[i] << endl;
        cout << endl;
    }

    cout << "Total Sequence Count: " << seq_count << endl;
    cout << "Length of the Alignment: " << ultimate_center.length() << endl;
    cout << endl;

    cout << "-------------------------------------------------------" << endl;
    cout << endl;

    cout << "Now let me write those stuff to a file to not forget!" <<endl;

    cout << endl;

    ofstream output(output_phy);
    output << "\t" << seq_count << "\t" << ultimate_center.length() << "\n";

    string finals[seq_count];

    int ind = 0;
    for(int i = 0; i < seq_count; i++){
        if(i == min_index){
            finals[i] = ultimate_center;
        } else {
            finals[i] = aligned_pairs[ind++];
        }
    }

    size_t line_count = (ultimate_center.length() / 50) + 1;

    size_t alignment_length = ultimate_center.length();
    string writing[line_count * seq_count];

    for(size_t i = 0; i < line_count * seq_count; i++){
        writing[i] = "";
    }

    for(int i = 0; i < seq_count; i++){
        size_t p = 0;
        int pos = i;

        while(seq_names[i].length() < 10){
            seq_names[i] += " ";
        }

        writing[pos] += seq_names[i];

        while(p < alignment_length){
            if(p > 0 && p % 50 == 0){
                pos = pos + seq_count;
                writing[pos] += "          ";
            }

            else if(p > 0 && p % 10 == 0)
                writing[pos] += " ";

            writing[pos] += finals[i].at(p);
            p++;
        }
    }

    for(size_t i = 0; i < line_count * seq_count; i++){
        output << writing[i] << endl;

        if( (i + 1) % seq_count == 0){
            output << endl;
        }
    }


    output.close();

    cout << "All done!" << endl;

    cout << endl;
    cout << "-------------------------------------------------------" << endl;

    cout << "BBBBBBBBBBBBBYYYYYYYYYYYYYYYYYYYYYEEEEEEEEEEEEEEEEEEEE" << endl;


    return 0;
}
