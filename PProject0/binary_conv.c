#include <stdio.h>
#include <math.h>

int i = 0;
int sum = 0;

void reverse(void)
{
    char c;
    if((c = getchar()) != '\n'){ 
        reverse();
        sum += (c -'0')*pow(2,i);
        i++;
    }
    return;
}

int main(void)
{
    printf("Enter a Binary number:\n");
    reverse();
    printf("The Decimal is: %d\n",sum);
    return ;
}