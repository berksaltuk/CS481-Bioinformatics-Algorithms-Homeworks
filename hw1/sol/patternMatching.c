#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>
#define _GNU_SOURCE
#define MIN(a,b) (((a)<(b))?(a):(b))

int nodeCount = 1;

char* substr(const char *src, int m, int n)
{

    int len = n - m;
    char *dest = (char*)malloc(sizeof(char) * (len + 1));

    for (int i = m; i < n && (*(src + i) != '\0'); i++)
    {
        *dest = *(src + i);
        dest++;
    }

    *dest = '\0';

    return dest - len;
}

//////////////////////////////////////////////////////////////////////////
////////////// Dynamic Array implementation - begin /////////////////////

struct intArray {

    int * arr;
    int cur_size;
    int actual_size;

};

void initIntArray(struct intArray *a, int size){
    a->arr = malloc(size * sizeof(int));
    a->cur_size = 0;
    a->actual_size = size;
}

void insertToArray(struct intArray * a, int val){
    if(a->cur_size == a->actual_size){
        a->actual_size *= 2;
        a->arr = realloc(a->arr, a->actual_size * sizeof(int));
    }
    a->arr[a->cur_size] = val;
    a->cur_size++;
}

void freeIntArray(struct intArray * a){
    free(a->arr);

    a->arr = NULL;
}

/////////////////////////////////////////////////////////////////////////
////////////// Keyword Tree implementation - begin /////////////////////

#define LETTER_COUNT 26

struct KeywordNode {
    struct KeywordNode *descendants[LETTER_COUNT];
    int nodeNo;
    bool last;
    struct intArray indexes;
};

struct KeywordNode *newKeywordNode(int label){
    struct KeywordNode *keywordNode = (struct KeywordNode *) malloc(sizeof(struct KeywordNode));
    for(int i = 0; i < LETTER_COUNT; i++){
        keywordNode->descendants[i] = NULL;
    }
    initIntArray(&keywordNode->indexes, 10);
    insertToArray(&keywordNode->indexes, label);
    return keywordNode;
};

void insertKeywordNode(struct KeywordNode *root, const char *text, int label){

    int len = strlen(text);

    struct KeywordNode *cur = root;
    for(int i = 0; i < len; i++){
        int letter_i = text[i] - 'a';

        if(cur->descendants[letter_i] == NULL) {
            cur->descendants[letter_i] = newKeywordNode(label);
            cur = cur->descendants[letter_i];
        }

        else {
            cur = cur->descendants[letter_i];

            insertToArray(&cur->indexes, label);
        }

    }
};

bool match(struct KeywordNode * root, const char *pattern, int index, FILE * f){

    int len = strlen(pattern);

    struct KeywordNode *cur = root;

    for(int i = 0; i < len; i++){
        int letter_i = pattern[i] - 'a';

        if(cur->descendants[letter_i] == NULL) {

            return false;
        }
        cur = cur->descendants[letter_i];
    }

    if(f){

        fprintf(f, "%d: ", index);
        for(int i = 0; i < cur->indexes.cur_size; i++){
            fprintf(f, "%d ", cur->indexes.arr[i]);
        }
        fprintf(f, "\n");

    } else {
        printf("%d: ", index);
        for(int i = 0; i < cur->indexes.cur_size; i++){

            printf("%d ", cur->indexes.arr[i]);
        }
        printf("\n");
    }

    return true;
}

void createAllSub(struct KeywordNode * root, const char * text){

    int textLen = strlen(text);

    for(int i = 0; i < textLen; i++){

        char * sub = substr(text, i, textLen);
        insertKeywordNode(root, sub, i + 1);

    }
}

void freeKeywordTree(struct KeywordNode * root){
    if(root){
        for(int i = 0; i < LETTER_COUNT; i++){
            freeKeywordTree(root->descendants[i]);
        }
        freeIntArray(&root->indexes);
        free(root);
    }
}

////////////// Keyword Tree implementation - end ////////////////////////
////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////
////////////// Suffix Tree implementation - begin //////////////////////

struct SuffixNode {
    struct SuffixNode * descendants[LETTER_COUNT + 1];
    int nodeCount;
    int suffixNo;
    char *suffix;
};

struct SuffixNode * newSuffixNode(char * suffix, int suffixNo, int nodeCount){
    struct SuffixNode * suffixNode = (struct SuffixNode * ) malloc(sizeof(struct SuffixNode));

    for(int i = 0; i < LETTER_COUNT + 1; i++){
        suffixNode->descendants[i] = NULL;
    }
    suffixNode->suffix = suffix;
    suffixNode->nodeCount = nodeCount;
    suffixNode->suffixNo = suffixNo;
    return suffixNode;
};


void insertSuffixNode(struct SuffixNode * root, char * given, int suffixNo){

    struct SuffixNode * cur = root;

    int givenLen = strlen(given);

    int suffixIndex = 0;
    int i = 0;
    while( i < givenLen ){
        if(cur->descendants[(given[i] - 'a')] == NULL){

            char * sub= substr(given, i, givenLen);
            char * subWD = malloc(strlen(sub)+1+1);
            strcpy(subWD, sub);
            strcat( subWD, "$");

            free(sub);

            cur->descendants[given[i] - 'a'] = newSuffixNode(subWD, suffixNo, nodeCount++);
            break;
        } else {

            cur = cur->descendants[(given[i] - 'a')];
            if( suffixIndex < strlen(cur->suffix) && cur->suffix[suffixIndex] == given[i]) {
                while( i < givenLen && suffixIndex < strlen(cur->suffix) && cur->suffix[suffixIndex] == given[i]){
                    i++;
                    suffixIndex++;
                }
            }
            if(i == givenLen && suffixIndex < strlen(cur->suffix)){
                char * com = substr(cur->suffix, 0, suffixIndex);

                char * rem = substr(cur->suffix, suffixIndex, strlen(cur->suffix));

                char * dol = "$";
                free(cur->suffix);

                cur->suffix = malloc(strlen(com) + 1);

                strcpy(cur->suffix, com);

                free(com);

                struct SuffixNode * newChild = newSuffixNode(rem, cur->suffixNo, nodeCount++);

                for(int i = 0; i < LETTER_COUNT + 1; i++){
                    if(cur->descendants[i]){
                        newChild->descendants[i] = cur->descendants[i];
                        cur->descendants[i] = NULL;
                    }
                }
                cur->descendants[rem[0] - 'a'] = newChild;

                struct SuffixNode * dollar = newSuffixNode(dol, suffixNo, nodeCount++);

                cur->descendants[LETTER_COUNT] = dollar;
                break;
            }
        }
    }

}

void addAllSuffixes(struct SuffixNode * root, char * text){

    int textLen = strlen(text);

    for(int i = 0; i < textLen; i++){
        char * sub = substr(text, i, textLen);
        insertSuffixNode(root, sub, i + 1);
    }
}


void findAllSuffixMatches(struct SuffixNode * cur){

    if(cur->suffix[strlen(cur->suffix) - 1] == '$'){
        printf("%d \n", cur->suffixNo);
    } else {
        for(int i = 0; i < LETTER_COUNT + 1; i++){
            if(cur->descendants[i])
                findAllSuffixMatches(cur->descendants[i]);
        }
    }

}

void matchSuffix(struct SuffixNode * root, char * key){

    struct SuffixNode* cur = root->descendants[key[0] - 'a'];

    int suffixIndex = 0;
    int i = 0;
    while(i < strlen(key)){
        if(!cur){
            return;
        }
        else {
            while(i < strlen(key) && suffixIndex < strlen(cur->suffix) && key[i] == cur->suffix[suffixIndex])
            {
                i++;
                suffixIndex++;
            }
            if(i == strlen(key))
            {
                findAllSuffixMatches(cur);
                break;
            }
            else if(suffixIndex == strlen(cur->suffix)){
                cur = cur->descendants[key[i] - 'a'];
                suffixIndex = 0;
            }
            else {
                return;
            }
        }
    }

    return;

}

void freeSuffixTree(struct SuffixNode * root){
    if(root){
        for(int i = 0; i < LETTER_COUNT + 1; i++){
            freeSuffixTree(root->descendants[i]);
        }

        free(root);
    }
}

void displayDOT(struct SuffixNode * root, FILE * f){

    bool hasDollar = false;
    for(int i = 0; i < strlen(root->suffix); i++){
        if(root->suffix[i] == '$'){
            hasDollar = true;
            break;
        }
    }

    if(hasDollar) {
        fprintf(f, "n%d [label=\"%d\" shape=box];\n", root->nodeCount, (root->suffixNo - 1));
    }
    else {
        fprintf(f, "n%d [label=\"\"];\n", root->nodeCount);
    }

    for(int i = 0; i < LETTER_COUNT + 1; i++){
        if(root->descendants[i]){
            fprintf(f, "n%d -> n%d [label=\"%s\"];\n", root->nodeCount, root->descendants[i]->nodeCount, root->descendants[i]->suffix);
        }
    }

    for(int i = 0; i < LETTER_COUNT + 1; i++){
        if(root->descendants[i]){
            displayDOT(root->descendants[i], f);
        }
    }

}

////////////// Suffix Tree implementation - end /////////////////////////
////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {

     if (argc != 7)
    {
        printf("You did not give all arguments!");
        exit(1);
    }

    if (strcmp(argv[1], "-p") != 0 || strcmp(argv[3], "-t") != 0 || strcmp(argv[5], "-o") != 0)
    {
        printf("Some of switches -p, -t or -o missing!");
        exit(2);
    }

    struct KeywordNode *root = newKeywordNode(0);
    struct SuffixNode * root2 = newSuffixNode("", 0, nodeCount++);

    FILE *fp = fopen(argv[4], "r");
    if(fp == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    char line[128];
    while(fgets(line, sizeof(line), fp) != NULL) {
        line[strcspn(line, "\n")] = 0;

        printf("Adding text %s to both trees...\n", line);
        createAllSub(root, line);
        addAllSuffixes(root2, line);
    }

    fclose(fp);

    printf("\n");

    printf("- - - KEYWORD TREE - - -\n");

    FILE *fp1 = fopen(argv[2], "r");
    if(fp1 == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    int i = 1;
    while(fgets(line, sizeof(line), fp1) != NULL) {
        line[strcspn(line, "\n")] = 0;
        printf("Matching: %s\n", line);
        if(!match(root, line, i++, NULL)){
            printf("No matches for pattern %d", i);
        }
    }

    fclose(fp1);

    printf("\n");
    printf("- - - SUFFIX TREE - - -\n");

    FILE *fp2 = fopen(argv[2], "r");
    if(fp2 == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    while(fgets(line, sizeof(line), fp2) != NULL) {
        line[strcspn(line, "\n")] = 0;
        printf("Matching: %s\n", line);
        matchSuffix(root2, line);
    }

    fclose(fp2);

    printf("\n");
    printf("Writing matching output to %s file...\n", argv[6]);

    FILE *f = fopen(argv[6], "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    FILE *fp3 = fopen(argv[2], "r");
    if(fp3 == NULL) {
        perror("Unable to open file!");
        exit(1);
    }

    i = 1;
    while(fgets(line, sizeof(line), fp3) != NULL) {
        line[strcspn(line, "\n")] = 0;
        if(!match(root, line, i++, f)){
            fprintf(f, "No matches for pattern %d", i);
        }
    }

    fclose(fp3);

    fclose(f);

    printf("Writing DOT file for Suffix Tree...\n");
    FILE *f2 = fopen("output.dot", "w");
    if (f2 == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    fprintf(f2, "diagraph BST {\n");

    struct SuffixNode * cur = root2;
    displayDOT(cur, f2);

    fprintf(f2, "}");
    fclose(f2);

    freeKeywordTree(root);
    freeSuffixTree(root2);

    return 0;
}
