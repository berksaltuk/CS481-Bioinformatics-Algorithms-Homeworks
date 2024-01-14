#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

struct MNode
{
    int row;
    int col;
    int val;

    struct MNode *prev;
    struct MNode *next;
};

struct MNode *head1 = NULL;
struct MNode *head2 = NULL;
struct MNode *head3 = NULL;
struct MNode *head4 = NULL;

void create(struct MNode **first, int n)
{
    struct MNode *cur;

    int total_nodes = 0;
    int row = 0;

    while (total_nodes <= n * n)
    {

        int rand_value = ((rand() % 1999) + (-999));

        if (rand_value != 0)
        {
            if ((*first) == NULL)
            {
                (*first) = (struct MNode *)malloc(sizeof(struct MNode));
                (*first)->next = (*first);
                (*first)->prev = (*first);

                (*first)->col = total_nodes % n;
                (*first)->row = row;
                (*first)->val = rand_value;
            }
            else
            {

                cur = (*first)->prev;

                struct MNode *mnode = (struct MNode *)malloc(sizeof(struct MNode));

                cur->next = mnode;
                mnode->prev = cur;
                mnode->next = (*first);
                (*first)->prev = mnode;

                mnode->col = total_nodes % n;
                mnode->row = row;

                mnode->val = rand_value;
            }
        }
        total_nodes++;
        if (total_nodes % n == 0)
        {
            row++;
        }
    }
}

void traverse(struct MNode *first, int n)
{
    struct MNode *cur = first;
    int total_nodes = 0;

    while (total_nodes < n * n)
    {
        int index = cur->row * n + cur->col;

        if (index > total_nodes)
        {
            printf("%6d", 0);
        }
        else
        {
            printf("%6d", cur->val);
            cur = cur->next;
        }
        total_nodes++;

        if (total_nodes % n == 0)
        {
            printf("\n");
        }
    }
}

void destroy(struct MNode **first)
{
    while (*first && (*first)->next != (*first))
    {
        struct MNode *temp = (*first)->next;
        temp->next->prev = (*first);
        (*first)->next = temp->next;

        free(temp);
    }

    free(*first);
    first = NULL;
}

void matrix_sum(struct MNode *first1, struct MNode *first2, int n)
{
    struct MNode *cur1 = first1;
    struct MNode *cur2 = first2;

    int total_nodes = 0;
    int i1, i2;

    while (total_nodes < n * n)
    {
        i1 = cur1->row * n + cur1->col;
        i2 = cur2->row * n + cur2->col;

        if (total_nodes < i1 && total_nodes < i2)
        {
            printf("%6d", 0);
        }
        else if (total_nodes == i1 && total_nodes == i2)
        {
            printf("%6d", cur1->val + cur2->val);
            cur1 = cur1->next;
            cur2 = cur2->next;
        }
        else if (total_nodes == i1)
        {
            printf("%6d", cur1->val);
            cur1 = cur1->next;
        }
        else
        {
            printf("%6d", cur2->val);
            cur2 = cur2->next;
        }

        total_nodes++;

        if (total_nodes % n == 0)
        {
            printf("\n");
        }
    }
}

void scalar_mult(struct MNode *first, int s)
{
    struct MNode *cur = first;
    while (cur->next != first)
    {
        cur->val = cur->val * s;
        cur = cur->next;
    }
}

void transpose(struct MNode *first, int n)
{
    struct MNode *cur = first;

    int found = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            while (cur->next != first)
            {

                if (cur->col == i && cur->row == j)
                {
                    printf("%6d", cur->val);
                    found = 1;
                    break;
                }
                cur = cur->next;
            }

            if (!found)
            {
                printf("%6d", 0);
            }
            cur = first;
            found = 0;
        }

        printf("\n");
    }
}

int main(int argc, char *argv[])
{

    int n = 0;
    int s = 0;
    if (argc != 5)
    {
        printf("You did not give 5 arguments!");
        exit(1);
    }

    if (strcmp(argv[1], "-n") != 0 || strcmp(argv[3], "-s") != 0)
    {
        printf("Either switch -n or -s is missing!");
        exit(2);
    }

    if (atoi(argv[2]) < 1)
    {
        printf("Dimension must be at least 1!");
        exit(3);
    }

    n = atoi(argv[2]);
    s = atoi(argv[4]);

    srand(time(NULL)); // seeding

    create(&head1, n); // creating n*n matrix
    create(&head2, n);

    printf("A:\n");
    traverse(head1, n);

    printf("B: \n");
    traverse(head2, n);

    printf("C: \n");
    matrix_sum(head1, head2, n);

    destroy(&head1);
    destroy(&head2);

    printf("\n");

    create(&head3, n);

    printf("s: %2d\n", s);

    printf("A: \n");
    traverse(head3, n);

    scalar_mult(head3, s);
    printf("C: \n");
    traverse(head3, n);

    destroy(&head3);

    create(&head4, n);
    printf("A: \n");
    traverse(head4, n);
    printf("A.T: \n");
    transpose(head4, n);
    destroy(&head4);
    return 0;
}