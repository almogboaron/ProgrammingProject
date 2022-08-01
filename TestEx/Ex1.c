#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
/**/
#define f(x,y) (x+=y;y=x-y;x-=y)
/**/
typedef struct node {
    int val;
    struct node* left;
    struct node* right;
    } Node;
Node* createNode(int);
Node* InsertNode(Node* root, int val){
    if(root==NULL){
        return createNode(val);
    }
    else if ( val <= root->val){
        root->left = InsertNode(root->left,val);
        return root;
    }else if( root->val < val){
        root->right = InsertNode(root->right,val);
        return root;
    }
}
Node* createNode(int val){
    Node* node = (Node*) malloc(sizeof(Node));
    node->val = val;
    node->left = NULL;
    node->right = NULL;
    return node;
}
/**/
typedef struct Stack{
    struct stackNode* first;
}Stack
typedef struct stackNode{
    int val ;
    stackNode* prev;
}stackNode;
typedef struct product{
    int id;
    double price;
}product;

Stack* init(int* arr, int n){
    int i;
    Stack* stack = (Stack*) malloc(sizeof(Stack));
    stackNode* prev = NULL;
    for(i=n-1;i>=0;i--){
        stackNode* node = (stackNode*) malloc(sizeof(stackNode));
        node->val = arr[i];
        node->prev = prev;
        prev = node;
    }
    stack->first = node;
    return stack;
}
int pop(Stack* stack){
    int val;
    stackNode* node = stack->first;
    stack->first = stack->first->prev;
    val = node->val;
    free(node);
    return val;
}
/**/
void sequence_set(unsigned int*flags,unsigned int i,unsigned int j){
    unsigned int size,number;
    size = sizeof(unsigned int)*8;
    for(;i<=j;i++){
        number = 1<<(i%size);
        arr[i/size]|= arr[i/size]|number;
    }
}
/**/
void readMat(*char path,int n){
    int *data,temp,**matrix,i,j,cols,rows,sum;
    FILE *ptf = fopen(path,r);
    data = (int*) malloc(sizeof(int)n*n);
    temp = fread(data,sizeof(int),n*n,ptf);
    asser(temp==n*n);

    matrix = (int**) malloc(sizeof(int*)*n);
    assert(matrix!=NULL);
    for(i=0;i<n;i++){
        matrix[i] = data +i*n;
    }
    sum=0;
    for(i=0;i<n;i++){
        cols = 0;
        rows = 0;
        for(j=0;j<n;j++){
            cols += matrix[j][i];
            rows += matrix[i][j];
        }
        sum+= rows>cols;
    }
    print(sum);
    free(data);
    free(matrix);

}

int my_compare(cont *void product1,const *void product2){
    if(p1->id == product2->id){return 0;}
    return product1->id > product2->id ? -1:1;
}
int main(int argc , char* argv[]){
    /*
    Node* root = createNode(4);
    printf("Hello");
    root = InsertNode(root,5);
    root = InsertNode(root,3);
    root = InsertNode(root,3);
    printf("%d %d %d %d",root->val,root->left->val,root->left->left->val,root->right->val);*/

}