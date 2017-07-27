#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_LINE_LENGTH 200 //Max no. of character of a line in the data file

int DIM_X=65;
int DIM_Y=65;
int DIM_Z=65;


int ReadLineSS(FILE*,double*,char*);  
int Get_Basic_Parameters(char in_file[]);
int Generate_vtk(char in_file[],char out_file[]);
int GetDIM(char in_file[]);  

int REPORT=0;
int SORTED=0;
main(){
      char s1[100],s2[100];
      int i,j,k;
      printf("Input filename?>");
      scanf("%s",s1);
      getchar();
      printf("Output filename?>");
      scanf("%s",s2);
      getchar();
      GetDIM(s1);
      Generate_vtk(s1,s2);
      printf("\nFINISHED\n");         
}

int ReadLineSS(FILE* fp_in,double data[],char line[]){
    /*
    Return value:
           0:EOF
           1:Data
           -1:Data without enough column
           2:time data
           -100:Empty line
    */
    int i,j,k,result;
    double f;
    char s[MAX_LINE_LENGTH];
    char s1[MAX_LINE_LENGTH];
    char s2[MAX_LINE_LENGTH];
    char s3[MAX_LINE_LENGTH];
    char s4[MAX_LINE_LENGTH];
    char s5[MAX_LINE_LENGTH];
    char s6[MAX_LINE_LENGTH];
    //Get a line of data using fgets
    if(fgets(line,MAX_LINE_LENGTH,fp_in)==NULL){
      return 0;                             
    }
    sscanf(line,"%s %s %s %s %s %s",s1,s2,s3,s4,s5,s6);
    
    if(s1[0]=='t'){
      //It is a line indicating the time
      f=atof(s2);
      data[0]=f;
      return -1;              
    }else if(s1[0]=='#'){
      //It is a comment line
      return -2;              
    }else if(line[0]=='\n' || s1[0]=='\n'){
      //It is an empty line
      return -100;    
    }else{
      //Probably it is a data line
      f=atof(s1);
      data[0]=f;
      f=atof(s2);
      data[1]=f;
      f=atof(s3);
      data[2]=f;
      f=atof(s4);
      data[3]=f;
      f=atof(s5);
      data[4]=f;
      f=atof(s6);
      data[5]=f;
    }
    return 6;
}

int GetDIM(char in_file[]){   
         
    //The data file need to be sorted

    
    FILE* fp_in;
    FILE* fp_out;
    int i;
    double nowxyz[3],firstxyz[6],lastxyz[6];
    
    int linecnt;
    int xcnt,ycnt,zcnt;
    double data[6];
    char filename[100];
    char line[MAX_LINE_LENGTH];
    int return_value;
    int flagx,flagy,flagz;
    
    //OPEN FILE
    fp_in=fopen(in_file,"r");
    if(fp_in==NULL){
      printf("INPUT FILE NOT EXIST\n");
      getchar();
      return 0;     
    }
    //READ FILE
    linecnt=0;
    xcnt=0;
    ycnt=0;
    zcnt=0;
    flagx=0;
    flagy=0;
    flagz=0;
    while(1){
    return_value=ReadLineSS(fp_in,data,line);
      //PROCESS FILE
      if(return_value>0){
        if(linecnt==0){
         nowxyz[0]=data[0];
         nowxyz[1]=data[1];
         nowxyz[2]=data[2];
         firstxyz[0]=data[0];
         firstxyz[1]=data[1];
         firstxyz[2]=data[2];
         xcnt++;
         ycnt++;
         zcnt++;
         linecnt++;
        }else{
           
           if(nowxyz[0]!=data[0] && !flagx){
             if(data[0]==firstxyz[0]){
               flagx=1;            
             }else{
               xcnt++;                  
               nowxyz[0]=data[0];
             }     
           }
           if(nowxyz[1]!=data[1] && !flagy){
             if(data[1]==firstxyz[1]){
               flagy=1;            
             }else{
               ycnt++;                  
               nowxyz[1]=data[1];
             }     
           }
           if(nowxyz[2]!=data[2] && !flagz){
             if(data[2]==firstxyz[2]){
               flagz=1;            
             }else{
               zcnt++;                  
               nowxyz[2]=data[2];
             }     
           }
           linecnt++;
           if(flagx && flagy && flagz){
             break;         
           }
        }
      }else if(return_value==0){
        break;      
      }
    }
    fclose(fp_in);
    DIM_X=xcnt;
    DIM_Y=ycnt;
    DIM_Z=zcnt;
    return 1; 
}


int Generate_vtk(char in_file[],char out_file[]){
    //
    
   int TOTAL_NO_OF_DATA; 
   TOTAL_NO_OF_DATA=(DIM_X)*(DIM_Y)*(DIM_Z); 
   FILE* fp_in;
   FILE* fp_in_rt;
   FILE* fp_out;
   FILE* fp_out2;
   FILE* fp_out3;
   FILE* fp_out4;
   //
    int i;
    int xcnt,ycnt,zcnt;
    int spacecnt,linecnt;
    double data[6];
    double data_rt[6];
    char c;
    char datatype[10]="double";
    char temp_file[100]="vtk1.tmp";
    char line[MAX_LINE_LENGTH];
    int return_value,return_value2;
    double firstxyz[6],lastxyz[6];
    //OPEN FILE
    fp_in=fopen(in_file,"r");
    if(fp_in==NULL){
      printf("INPUT FILE NOT EXIST\n");
      return 0;     
    }
    fp_out=fopen(out_file,"w");
    fp_out2=fopen(temp_file,"w");
    //
    if(fp_out==NULL){
      printf("OPEN OUTPUT FILE FAIL\n");

      return 0;     
    }
    //READ FILE
    spacecnt=1;
    linecnt=0;
    //WRITE PREFIX
    fprintf(fp_out,"# vtk DataFile Version 3.0\n");
    fprintf(fp_out,"vtk output\n");
    fprintf(fp_out,"ASCII\n");
    
    fprintf(fp_out,"DATASET STRUCTURED_GRID\n");
    fprintf(fp_out,"DIMENSIONS %d %d %d\n",DIM_Z,DIM_Y,DIM_X);      
    
    fprintf(fp_out,"POINTS %d float\n",TOTAL_NO_OF_DATA);
    //
    fprintf(fp_out2,"POINT_DATA %d\n",TOTAL_NO_OF_DATA);
    fprintf(fp_out2,"VECTORS curr_density %s\n",datatype);
    //WRITE DATA
    while(1){
      return_value=ReadLineSS(fp_in,data,line);
      if(return_value==0){
        break;                    
      }else if(return_value>0){
            fprintf(fp_out,"%.2f\t%.2f\t%.2f\n",data[0],data[1],data[2]); 
            fprintf(fp_out2,"%.5E\t%.5E\t%.5E\n",data[3],data[4],data[5]);  
        linecnt++;
      }
    }
    fclose(fp_in);
    fclose(fp_out);
    fclose(fp_out2);
    //Combine the files
    
    printf("--NOW COMBINING FILES--\n");
    fp_out=fopen(out_file,"a");
    fp_out2=fopen(temp_file,"r");
    //
    //
    while((c=fgetc(fp_out2))!=EOF){
      fputc(c, fp_out);
    }
    //
    fclose(fp_out);
    fclose(fp_out2);
    
    return 1;
}


//QUICK SORT
void swap(double *a, double *b)
{
  double t=*a; *a=*b; *b=t;
}
void quick_sort(double arr[], int beg, int end)
{
  if (end > beg + 1)
  {
    double piv = arr[beg]; 
    int l = beg + 1, r = end;
    while (l < r)
    {
      if (arr[l] <= piv)
        l++;
      else
        swap(&arr[l], &arr[--r]);
    }
    swap(&arr[--l], &arr[beg]);
    quick_sort(arr, beg, l);
    quick_sort(arr, r, end);
  }
}
