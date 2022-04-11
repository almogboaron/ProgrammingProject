#include <stdio.h>
#include <math.h>

int i = 0;
int sum = 0;
int a = 0;
int b = 0;
int num = 0;

void reverse(void)
{
    char c;
    convertBase
    if((c = getchar()) != '\n'){ 
        reverse();
        sum += (c -'0')*pow(2,i);
        i++;
    }
    return;
}

int main(void)
{
    //from base a
    printf("Please enter the numbers base:\n");
    scanf("%d",&a);
    if(!(2 <= a <= 16)) {
        printf("Invalid input base");
        exit();
    }
    
    //To base b
    printf("Please enter the desired base:\n");
    scanf("%d",&b);
    if(!(2 <= b <= 16)) {
        printf("Invalid input base");
        exit();
    }
    
    //Changing the number num:
    printf("Please enter a number in base %d\n",a);
    convertBase();


    printf("The Decimal is: %d\n",sum);
    return ;
}