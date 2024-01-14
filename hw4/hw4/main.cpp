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

double calculate_edit_distance(string seq1, string seq2, int match, int mismatch_penalty, int gap_penalty){

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

tuple<int, int> find_min_distance(vector<vector<double>> matrix, double &height){

    double distance = INT_MAX;

    for(size_t i = 0; i < matrix.size(); i++){
        for(size_t j = 0; j < i; j++){
            if(matrix[i][j] < distance)
            {
                distance = matrix[i][j];
            }
        }
    }
    height = distance/2;

    int i_min = -1, j_min = -1;
    for(size_t i = 0; i < matrix.size(); i++){
        for(size_t j = 0; j < i; j++){
            if(matrix[i][j] == distance)
            {
                i_min = i;
                j_min = j;
                return make_tuple(i_min, j_min);
            }
        }
    }

    return make_tuple(i_min, j_min);
}

double calculate_new_distance(vector<vector<double>> matrix, int i_min, int j_min, int i, vector<vector<string>> sequences){

    if(i < i_min && i < j_min){
        return (double)(matrix[i_min][i] * sequences[i_min].size() + matrix[j_min][i] * sequences[j_min].size())/
            (sequences[i_min].size() + sequences[j_min].size() );
    }
    else if(i < i_min && i > j_min){
        return (double)(matrix[i_min][i] * sequences[i_min].size() + matrix[i][j_min] * sequences[j_min].size())/
            (sequences[i_min].size() + sequences[j_min].size());
    }
    else if(i > i_min && i < j_min){
        return (double)(matrix[i][i_min] * sequences[i_min].size() + matrix[j_min][i] * sequences[j_min].size())/
            (sequences[i_min].size() + sequences[j_min].size());
    }
    else{
        return (double)(matrix[i][i_min] * sequences[i_min].size() + matrix[i][j_min] * sequences[j_min].size())/
            (sequences[i_min].size() + sequences[j_min].size());
    }

}

vector<vector<double>> update_distances(vector<vector<double>> matrix, int i_min, int j_min, vector<vector<string>> &sequences){

    int high = max(i_min, j_min);
    int low = min(i_min, j_min);

    vector<vector<double>> new_matrix = matrix;

    for(int i = 0; i < new_matrix.size(); i++){

        for(int j = 0; j < i; j++){

            if(j == low){
                new_matrix[i][j] = calculate_new_distance(matrix, i_min, j_min, i, sequences);
            }
            else if(i == low){
                new_matrix[i][j] = calculate_new_distance(matrix, i_min, j_min, j, sequences);

            } else {
                new_matrix[i][j] = matrix[i][j];

            }
        }
    }

    if (new_matrix.size() > high)
        new_matrix.erase(new_matrix.begin() + high);

    for (size_t i = 0; i < new_matrix.size(); i++)
    {
        if (new_matrix[i].size() > high)
            new_matrix[i].erase(new_matrix[i].begin() + high);
    }

    for(string s : sequences[high]){
        sequences[low].push_back(s);
    }

    sequences.erase(sequences.begin() + high);

    return new_matrix;
}

vector<vector<double>> neighbor_joining_update_distances(vector<vector<double>> matrix, vector<vector<string>> &sequences, string &output_str){

    double L = matrix.size();

    vector<double> r_vector;

    for ( int i = 0; i < L; i++ ) {
        double sum = 0;
        for ( int j = 0; j < L; j++ ){
            sum += matrix[i][j];
        }
        double r = (L > 2) ? ((double)(1 / (L - 2)) * sum) : 0;
        r_vector.push_back(r);
    }
    /*cout << "r vector" << endl;
    for(int i = 0; i < r_vector.size(); i++){
        cout << r_vector[i] << " ";
    }
    cout << endl;*/

    double distance = INT_MAX;

    int opt_i, opt_j = 0;

    for ( int i = 0; i < L - 1; i++ ) {

        for ( int j = i + 1; j < L; j++ ){

            double temp_distance = matrix[i][j] - (r_vector[i] + r_vector[j]);
            if(temp_distance < distance){
                distance = temp_distance;
                opt_i = i, opt_j = j;
            }
        }
    }

    double height_i = (double)(matrix[opt_i][opt_j] + r_vector[opt_i] - r_vector[opt_j]) / 2, height_j = (double)(matrix[opt_i][opt_j] - r_vector[opt_i] + r_vector[opt_j]) / 2;
    //cout << "heights " << height_i << " " << height_j << endl;
    if(L == 2){
        opt_i = 0, opt_j = 1;
        height_i = matrix[1][0], height_j = matrix[1][0];

        if(sequences[opt_i].size() == 1){
              output_str = "(" + output_str + "," + sequences[opt_i][0] + ":" + to_string(height_i) + ")";
        }
        else if(sequences[opt_j].size() == 1){
              output_str = "(" + output_str + "," + sequences[opt_j][0] + ":" + to_string(height_j) + ")";
        } else {
            output_str = "(" + output_str + ")";
        }
    }
    else if(sequences[opt_i].size() == 1 && sequences[opt_j].size() == 1){
        output_str += "(" + sequences[opt_i][0] + ":" + to_string(height_i) + "," + sequences[opt_j][0] + ":" + to_string(height_j) + ")";
    }
    else if(sequences[opt_i].size() > 1 && sequences[opt_j].size() == 1){
        output_str += ":" + to_string(height_i) + "," + sequences[opt_j][0] + ":" + to_string(height_j);
        output_str = "(" + output_str + ")";
    }
    else if(sequences[opt_i].size() == 1 && sequences[opt_j].size() > 1){
        output_str += ":" + to_string(height_j) + "," + sequences[opt_i][0] + ":" + to_string(height_i);
        output_str = "(" + output_str + ")";
    }

    if(sequences[opt_i].size() > 1 && sequences[opt_j].size() > 1){
        int found = 0;
        for(int i = output_str.length() - 1; i > 1; i--){
            if(output_str[i] == '(' && output_str[i-1] != ','){
                found = i;
                break;
            }
        }
        string left = output_str.substr(0,found);
        string right = output_str.substr(found);
        output_str = left + "," + right;
    }

    int high = max(opt_i, opt_j), low = min(opt_i, opt_j);

    vector<vector<double>> new_matrix = matrix;

    for(int i = 0; i < new_matrix.size(); i++){

        for(int j = 0; j < new_matrix.size(); j++){

            if(j == low){
                new_matrix[i][j] = new_matrix[j][i];
            }
            else if(i == low){
                new_matrix[i][j] = (double)(matrix[opt_i][j] + matrix[opt_j][j] - matrix[opt_i][opt_j]) / 2;
            }
            else {
                new_matrix[i][j] = matrix[i][j];

            }
        }
    }

    if (new_matrix.size() > high)
        new_matrix.erase(new_matrix.begin() + high);

    for (size_t i = 0; i < new_matrix.size(); i++)
    {
        if (new_matrix[i].size() > high)
            new_matrix[i].erase(new_matrix[i].begin() + high);
    }

    for(string s : sequences[high]){
        sequences[low].push_back(s);
    }

    sequences.erase(sequences.begin() + high);

    return new_matrix;
}

int main(int argc, char **argv){

    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    string input_fasta, output_tree, type;
    cmdl("-i") >> input_fasta;
    cmdl("-o") >> output_tree;
    cmdl("-t") >> type;

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

    vector<vector<double>> distance_matrix(seq_count, vector<double> (seq_count, 0));

    for(int i = 0; i < seq_count; i++){
        for( int j = 0; j < seq_count; j++){
            distance_matrix[i][j] = calculate_edit_distance(seq_arr[i], seq_arr[j], match, mismatch_penalty, gap_penalty);
        }
    }


    /*cout << "DISTANCE MATRIX FOR CHOOSING THE CENTER SEQUENCE: " << endl;
    cout << endl;
    for (int i = 0; i < seq_count; i++)
       for (int j = 0; j < seq_count; j++)
          cout << distance_matrix[i][j] << " \n"[j == seq_count-1];

    cout << endl;*/

    string output_str = "";

    vector<vector<string>> sequences(seq_count);
    for( int i = 0; i < seq_count; i++){
        sequences[i].push_back(seq_names[i]);
    }

    vector<double> tree_h;

    if( type.compare("upgma") == 0)
    {

        vector<vector<double>> adjusted = distance_matrix;

        int i_min, j_min;

        double h;
        for(int cur = 1; cur < seq_count; cur++){

            tie(i_min, j_min) = find_min_distance(adjusted, h);

            if(sequences[i_min].size() == 1 && sequences[j_min].size() == 1){
                output_str += "(" + sequences[i_min][0] + ":" + to_string(h) + "," + sequences[j_min][0] + ":" + to_string(h) + ")";
            }
            else if(sequences[i_min].size() > 1 && sequences[j_min].size() == 1){
                output_str += ":" + to_string(h - tree_h[tree_h.size()-1]) + "," + sequences[j_min][0] + ":" + to_string(h) + ")";
                output_str = "(" + output_str;
            }
            else if(sequences[i_min].size() == 1 && sequences[j_min].size() > 1){
                output_str += ":" + to_string(h - tree_h[tree_h.size()-1]) + "," + sequences[i_min][0] + ":" + to_string(h) + ")";
                output_str = "(" + output_str;
            }
            else if(sequences[i_min].size() > 1 && sequences[j_min].size() > 1){
                int found = 0;
                for(int i = output_str.length() - 1; i > 1; i--){

                    if(output_str[i] == '(' && output_str[i-1] != ','){
                        found = i;
                        break;
                    }
                }
                string left = output_str.substr(0,found);
                string right = output_str.substr(found);
                output_str = left + ":" + to_string(h - tree_h[tree_h.size()-2]) + "," + right;
                output_str += ":" + to_string(h - tree_h[tree_h.size()-1]) + ")";
                output_str = "(" + output_str;
            }

            vector<vector<double>> after_processing = update_distances(adjusted, i_min, j_min, sequences);
            adjusted = after_processing;


            tree_h.push_back(h);
        }

        //cout << output_str << endl;
    } else if (type.compare("nj") == 0) {

        /*vector<vector<double>> mat {
            {0,4,5,10},
            {4,0,7,12},
            {5,7,0,9},
            {10,12,9,0}
        };
        int seq_count2 = 4;
        string seq_names2[4] = {"A", "B", "C", "D"};*/
        vector<vector<string>> sequences(seq_count);
        for( int i = 0; i < seq_count; i++){
            sequences[i].push_back(seq_names[i]);
        }

        /*for (int i = 0; i < distance_matrix.size(); i++)
               for (int j = 0; j < distance_matrix.size(); j++)
                  cout << distance_matrix[i][j] << " \n"[j == distance_matrix.size() - 1];*/

        for(int cur = 1; cur < seq_count; cur++){
            vector<vector<double>> temp = neighbor_joining_update_distances(distance_matrix, sequences, output_str);
            distance_matrix = temp;
        }

    } else {

        cout << "Unknown method!" << endl;
        exit(1);
    }

    ofstream output(output_tree);
    output << output_str;
    output.close();

    cout << "Please check \"" << output_tree << "\" for output." << endl;
    return 0;
}
