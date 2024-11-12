//
//  main.c
//  myKING
//
//  Created by Meng Huang on 6/18/21.
//  Edited by Meng Huang on 10/2/24.
//  Copyright © 2021 Meng Huang. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


void *getFileName(char *pathname){
    char *q;
//    char ch = "\\"; //for windows pathway
//    char ch = "/"; //for linux pathway
    q = strrchr(pathname, '/') + 1;
    return q;
}

int remove_temp(char *fileName_input1, char *fileName_input2, char *path_input2){
    char fileName_temp[1000] = "";
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".bcf");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".bcf.csi");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".vcf");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".log");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".bim");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".bed");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".fam");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".pos");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".pre");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".raw");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".plink2_std");
    remove(fileName_temp);
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,"king_chi.log");
    remove(fileName_temp);
    return 0;
}
       
int compress_bit(char *input1, char *input2, char *output){
    /*compress numeric format to binary files*/
    int i=0, j=0, k=0, n_row=0, temp1=0,temp2=0,temp3=0;
    int n_marker=0, row=0;
    FILE *fp_g=NULL, *fp_a1=NULL, *fp_a0=NULL, *fp_bim=NULL, *fp_fam=NULL, *fp_r=NULL;
    char fileName_temp[1000] = "";
    strcpy (fileName_temp,output);
    strcat (fileName_temp,".king.txt");
    if((fp_r=fopen(fileName_temp,"a"))==NULL){
        printf("The .king.txt file isn't exist!\n");
        return -1;
    }
    
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".pre");
    if((fp_a1=fopen(fileName_temp,"wb+"))==NULL){
        printf("The file pre can't be created!\n");
        return -1;
    }
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".pos");
    if((fp_a0=fopen(fileName_temp,"wb+"))==NULL){
        printf("The file pos can't be created!\n");
        return -1;
    }
//    strcpy (fileName_temp,output);
//    strcat (fileName_temp,".list");
//    if((fp_m=fopen(fileName_temp,"w"))==NULL){
//        printf("The file SNP list can't be created!\n");
//        return -1;
//    }
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".raw");
    if((fp_g=fopen(fileName_temp,"r"))==NULL){
        printf("The genotype file doesn't exist!\n");
        return -1;
    }
    /*get m_marker and row*/
    
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".bim");
    if((fp_bim=fopen(fileName_temp,"r"))==NULL){
        printf("The .bim file doesn't exist!\n");
        return -1;
    }
    
    struct gm{
        char chr[100];
        char rs[100];
        char gen[100];
        char pos[100];
        char ref[5];
        char alt[5];
    };
    struct gm a[2]={0};
    while (fscanf(fp_bim,"%s%s%s%s%s%s", a[0].chr, a[0].rs, a[0].gen, a[0].pos, a[0].ref, a[0].alt)>0) {
        n_marker++;
    }
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".fam");
    
    if((fp_fam=fopen(fileName_temp,"r"))==NULL){
        printf("The .fam file doesn't exist!\n");
        return -1;
    }
    
    struct fm{
        char chr[100];
        char rs[100];
        char gen[100];
        char pos[100];
        char ref[100];
        char alt[100];
    };
    struct fm b[2]={0};
    while (fscanf(fp_fam,"%s%s%s%s%s%s", b[0].chr, b[0].rs, b[0].gen, b[0].pos, b[0].ref, b[0].alt)>0) {
        fprintf(fp_r, "%s\t", b[0].rs);
        row++;
    }
    
       
//    printf("number of individual is %d\nnumber of marker is %d\n",row, n_marker);
    if(row<1){
        printf("There are 0 individual!\n");
        return -1;
    }
    if(n_marker<1000){
        printf("There are less than 1000 SNPs!\n");
        return -1;
    }
    int n_residue=0;
    if ((n_marker/64 * 64) != n_marker) {
        n_row = (n_marker/64) + 1;
    }
    else {
        n_row = n_marker/64;
    }
    n_residue = 64*(n_row-1);
    
    unsigned long long int *a1 = (unsigned long long int*)calloc(n_row*3, sizeof(unsigned long long int));
    unsigned long long int *a0 = (unsigned long long int*)calloc(n_row*3, sizeof(unsigned long long int));
    char *buf = (char*)calloc(row*n_marker, sizeof(char));
    char *temp = (char*)calloc(100, sizeof(char));
    fwrite(&n_row, sizeof(int), 1, fp_a0);
    fwrite(&n_marker, sizeof(int), 1, fp_a0);
        
    for (i = 0; i < n_marker+6; i++) {
        fscanf(fp_g,"%s" ,temp);
//        fprintf(fp_m, "%s\n", temp);
//        printf("%s\n",temp);
    }
    for (i = 0; i < n_row*2; i++) {
        a1[i] = 0;
        a0[i] = 0;
    }
    for (k = 0; k < 2; k++) {
        for (i = 0; i < 6; i++) {
            fscanf(fp_g,"%s" ,temp);
//            printf("%s\n",temp);
        }
        temp1 = n_row * k;
        temp2 = n_marker * k;
        for (j=0; j < n_row-1; j++) {
            temp3 = 64*j;
            for (i=0; i<64; i++) {
                fscanf(fp_g,"%s" ,temp);
                if (temp[0]=='0') {
                    a1[temp1+j] <<= 1;
                    a0[temp1+j] <<= 1;
                }
                else if (temp[0]=='1') {
                    a1[temp1+j] <<= 1;
                    a0[temp1+j] <<= 1;
                    a0[temp1+j] = a0[temp1+j] | 1;
                }
                else if (temp[0]=='2') {
                    a1[temp1+j] <<= 1;
                    a0[temp1+j] <<= 1;
                    a1[temp1+j] = a1[temp1+j] | 1;
                }
                else if (temp[0]=='N') {
                a1[temp1+j] <<= 1;
                a0[temp1+j] <<= 1;
                a0[temp1+j] = a0[temp1+j] | 1;
                a1[temp1+j] = a1[temp1+j] | 1;
//                    printf("%s\n",temp);
                }
            }
        }
        for (i = temp2+n_residue; i < temp2+n_marker; i++) {
            fscanf(fp_g,"%s" ,temp);
            if (temp[0]=='0') {
                a1[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] <<= 1;
            }
            else if (temp[0]=='1') {
                a1[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] = a0[temp1+n_row-1] | 1;
            }
            else if (temp[0]=='2') {
                a1[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] <<= 1;
                a1[temp1+n_row-1] = a1[temp1+n_row-1] | 1;
            }
            else if (temp[0]=='N') {
            a1[temp1+n_row-1] <<= 1;
            a0[temp1+n_row-1] <<= 1;
            a0[temp1+n_row-1] = a0[temp1+n_row-1] | 1;
            a1[temp1+n_row-1] = a1[temp1+n_row-1] | 1;
            }
        }
        
        
    }
    fwrite(a1, sizeof(unsigned long long int), n_row*2, fp_a1);
    fwrite(a0, sizeof(unsigned long long int), n_row*2, fp_a0);
    free(a1);
    free(a0);
    free(buf);
    free(temp);
    fclose(fp_g);
    fclose(fp_bim);
    fclose(fp_fam);
    fclose(fp_r);
    fclose(fp_a1);
    fclose(fp_a0);
    return 0;
}

int compress_bit2(char *input1, char *path_input2, char *output){
    /*compress numeric format to binary files*/
    int i=0, j=0, k=0, n_row=0, temp1=0,temp2=0,temp3=0;
    int n_marker=0, row=0;
    FILE *fp_g=NULL, *fp_a1=NULL, *fp_a0=NULL, *fp_bim=NULL, *fp_fam=NULL, *fp_r=NULL;
    char fileName_temp[1000] = "";
    strcpy (fileName_temp,output);
    strcat (fileName_temp,".king.txt");
    if((fp_r=fopen(fileName_temp,"a"))==NULL){
        printf("The .king.txt file isn't exist!\n");
        return -1;
    }
    
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".pre");
    if((fp_a1=fopen(fileName_temp,"wb+"))==NULL){
        printf("The file pre can't be created!\n");
        return -1;
    }
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".pos");
    if((fp_a0=fopen(fileName_temp,"wb+"))==NULL){
        printf("The file pos can't be created!\n");
        return -1;
    }

    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".raw");
    if((fp_g=fopen(fileName_temp,"r"))==NULL){
        printf("The genotype file doesn't exist!\n");
        return -1;
    }
    /*get m_marker and row*/
    
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".bim");
    if((fp_bim=fopen(fileName_temp,"r"))==NULL){
        printf("The .bim file doesn't exist!\n");
        return -1;
    }
    
    struct gm{
        char chr[100];
        char rs[100];
        char gen[100];
        char pos[100];
        char ref[5];
        char alt[5];
    };
    struct gm a[2]={0};
    while (fscanf(fp_bim,"%s%s%s%s%s%s", a[0].chr, a[0].rs, a[0].gen, a[0].pos, a[0].ref, a[0].alt)>0) {
        n_marker++;
    }
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".fam");
    
    if((fp_fam=fopen(fileName_temp,"r"))==NULL){
        printf("The .fam file doesn't exist!\n");
        return -1;
    }
    
    struct fm{
        char chr[100];
        char rs[100];
        char gen[100];
        char pos[100];
        char ref[100];
        char alt[100];
    };
    struct fm b[2]={0};
    while (fscanf(fp_fam,"%s%s%s%s%s%s", b[0].chr, b[0].rs, b[0].gen, b[0].pos, b[0].ref, b[0].alt)>0) {
        fprintf(fp_r, "%s\t", b[0].rs);
        row++;
    }
    
       
//    printf("number of individual is %d\nnumber of marker is %d\n",row, n_marker);
    if(row<1){
        printf("There are 0 individual!\n");
        return -1;
    }
    if(n_marker<1000){
        printf("There are less than 1000 SNPs!\n");
        return -1;
    }
    int n_residue=0;
    if ((n_marker/64 * 64) != n_marker) {
        n_row = (n_marker/64) + 1;
    }
    else {
        n_row = n_marker/64;
    }
    n_residue = 64*(n_row-1);
    
    unsigned long long int *a1 = (unsigned long long int*)calloc(n_row*3, sizeof(unsigned long long int));
    unsigned long long int *a0 = (unsigned long long int*)calloc(n_row*3, sizeof(unsigned long long int));
    char *buf = (char*)calloc(row*n_marker, sizeof(char));
    char *temp = (char*)calloc(100, sizeof(char));
    fwrite(&n_row, sizeof(int), 1, fp_a0);
    fwrite(&n_marker, sizeof(int), 1, fp_a0);
        
    for (i = 0; i < n_marker+6; i++) {
        fscanf(fp_g,"%s" ,temp);
//        fprintf(fp_m, "%s\n", temp);
//        printf("%s\n",temp);
    }
    for (i = 0; i < n_row*2; i++) {
        a1[i] = 0;
        a0[i] = 0;
    }
    for (k = 0; k < 2; k++) {
        for (i = 0; i < 6; i++) {
            fscanf(fp_g,"%s" ,temp);
//            printf("%s\n",temp);
        }
        temp1 = n_row * k;
        temp2 = n_marker * k;
        for (j=0; j < n_row-1; j++) {
            temp3 = 64*j;
            for (i=0; i<64; i++) {
                fscanf(fp_g,"%s" ,temp);
                if (temp[0]=='0') {
                    a1[temp1+j] <<= 1;
                    a0[temp1+j] <<= 1;
                }
                else if (temp[0]=='1') {
                    a1[temp1+j] <<= 1;
                    a0[temp1+j] <<= 1;
                    a0[temp1+j] = a0[temp1+j] | 1;
                }
                else if (temp[0]=='2') {
                    a1[temp1+j] <<= 1;
                    a0[temp1+j] <<= 1;
                    a1[temp1+j] = a1[temp1+j] | 1;
                }
                else if (temp[0]=='N') {
                a1[temp1+j] <<= 1;
                a0[temp1+j] <<= 1;
                a0[temp1+j] = a0[temp1+j] | 1;
                a1[temp1+j] = a1[temp1+j] | 1;
//                    printf("%s\n",temp);
                }
            }
        }
        for (i = temp2+n_residue; i < temp2+n_marker; i++) {
            fscanf(fp_g,"%s" ,temp);
            if (temp[0]=='0') {
                a1[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] <<= 1;
            }
            else if (temp[0]=='1') {
                a1[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] = a0[temp1+n_row-1] | 1;
            }
            else if (temp[0]=='2') {
                a1[temp1+n_row-1] <<= 1;
                a0[temp1+n_row-1] <<= 1;
                a1[temp1+n_row-1] = a1[temp1+n_row-1] | 1;
            }
            else if (temp[0]=='N') {
            a1[temp1+n_row-1] <<= 1;
            a0[temp1+n_row-1] <<= 1;
            a0[temp1+n_row-1] = a0[temp1+n_row-1] | 1;
            a1[temp1+n_row-1] = a1[temp1+n_row-1] | 1;
            }
        }
        
        
    }
    fwrite(a1, sizeof(unsigned long long int), n_row*2, fp_a1);
    fwrite(a0, sizeof(unsigned long long int), n_row*2, fp_a0);
    free(a1);
    free(a0);
    free(buf);
    free(temp);
    fclose(fp_g);
    fclose(fp_bim);
    fclose(fp_fam);
    fclose(fp_r);
    fclose(fp_a1);
    fclose(fp_a0);
    return 0;
}

int popCount(unsigned long long int i){
    i = i - ((i >> 1) & 0x5555555555555555);
    i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
    i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
    i = i + (i >> 8);
    i = i + (i >> 16);
    i = i + (i >> 32);
    return (int)i & 0x7f;
}

int get_n1234(unsigned long long int *a1, unsigned long long int *a0, unsigned long long int *w1, unsigned long long int *w0, const int n_row, double *xw){
    int j=0,temp=0;
    for (j = 0; j < n_row; j++) {
        xw[0] += popCount(a0[j] & w0[j]);
        temp += (popCount(a1[j] & w0[j]) + popCount(a0[j] & w1[j]));
        xw[1] += popCount(a1[j] ^ w1[j]);
        xw[2] += popCount(a0[j]);
        xw[3] += popCount(w0[j]);
    }
    xw[1] = xw[1] -temp;
    
    return 0;
}

int kin_robust(double *xw){
    int n5 = xw[2]<=xw[3]?xw[2]:xw[3];
    xw[4] = (xw[0] - xw[1]*2)/(xw[2] + xw[3]);
    xw[5] = (xw[0] - 2*xw[1])/(xw[2] + xw[3]);//kin within family
    xw[6] = 0.5 - (0.25*(xw[2] + xw[3] - xw[0]*2) + xw[1])/n5;//kin0 cross family
    
    
    return 0;
}

int file_check(char *input1, char *input2, char *path_input2, char *fileName_output){
        int n_marker=0;
        FILE *fp_bim1=NULL, *fp_bim2=NULL, *fp_r=NULL;
        char fileName_temp[100] = "";
    char fileName_temp1[100] = "";
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(fileName_output));
    strcat (fileName_temp,".snp_list.txt");
    if((fp_r=fopen(fileName_temp,"w"))==NULL){
        printf("The results file pos can't be created!\n");
        return -1;
    }
        
        /*check marker number before input KING*/
        
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
        strcat (fileName_temp,".bim");
        if((fp_bim1=fopen(fileName_temp,"r"))==NULL){
            printf("The SNP list of file1 doesn't exist!\n");
            return -1;
        }
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input2));
        strcat (fileName_temp,".bim");
        if((fp_bim2=fopen(fileName_temp,"r"))==NULL){
            printf("The SNP list of file2 doesn't exist!\n");
            return -1;
        }
    
        struct gm{
            char chr[100];
            char rs[100];
            char gen[100];
            char pos[100];
            char ref[5];
            char alt[5];
        };
        struct gm a[2]={0};
        struct gm b[2]={0};
        while (fscanf(fp_bim1,"%s%s%s%s%s%s", a[0].chr, a[0].rs, a[0].gen, a[0].pos, a[0].ref, a[0].alt)>0 && fscanf(fp_bim2,"%s%s%s%s%s%s", b[0].chr, b[0].rs, b[0].gen, b[0].pos, b[0].ref, b[0].alt)>0) {
            if(strncmp(a[0].rs, "chr", 3)!=0){
                strcpy (fileName_temp,"chr");
                strcat (fileName_temp,a[0].rs);
            }
            if(strncmp(b[0].rs, "chr", 3)!=0){
                strcpy (fileName_temp1,"chr");
                strcat (fileName_temp1,b[0].rs);
            }
            if(strcmp(fileName_temp, fileName_temp1)==0){
                fprintf(fp_r, "%s\n", a[0].rs);
                n_marker++;
            }
        }
               
    fclose(fp_bim1);
    fclose(fp_bim2);
    fclose(fp_r);
    return n_marker;
}

int KING(char *input1, char *input2, char *path_input2, char *output){
    FILE *fp_a11=NULL, *fp_a01=NULL, *fp_a12=NULL, *fp_a02=NULL, *fp_r=NULL;
    int m=2, i=0, j=0, n_row1=0, n_row2=0, n_row=0, n_marker=0;
    unsigned long long int temp1=0, temp2=0, temp3=0, temp4=0, temp5=0, temp6=0;
    char fileName_temp[100] = "";
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".pre");
    if((fp_a11=fopen(fileName_temp,"rb"))==NULL){
        printf("The file .pre can't be opened!\n");
        return -1;
    }
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input1));
    strcat (fileName_temp,".pos");
    if((fp_a01=fopen(fileName_temp,"rb"))==NULL){
        printf("The file .pos can't be opened!\n");
        return -1;
    }
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".pre");
    if((fp_a12=fopen(fileName_temp,"rb"))==NULL){
        printf("The file .pre can't be opened!\n");
        return -1;
    }
    strcpy (fileName_temp,path_input2);
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".pos");
    if((fp_a02=fopen(fileName_temp,"rb"))==NULL){
        printf("The file .pos can't be opened!\n");
        return -1;
    }
    strcpy (fileName_temp,output);
    strcat (fileName_temp,".king.txt");
    if((fp_r=fopen(fileName_temp,"a"))==NULL){
        printf("The results file pos can't be created!\n");
        return -1;
    }
    fread(&n_row1, sizeof(int), 1,fp_a01);
    fread(&n_row2, sizeof(int), 1,fp_a02);
    if (n_row1 != n_row2){
        printf("The SNP number of two input compressed files are different!\n");
        return -1;
    }
    n_row=n_row1;
    unsigned long long int *a1 = (unsigned long long int*)calloc(n_row1*m, sizeof(unsigned long long int));
    unsigned long long int *a0 = (unsigned long long int*)calloc(n_row1*m, sizeof(unsigned long long int));
    unsigned long long int *w1 = (unsigned long long int*)calloc(n_row2*m, sizeof(unsigned long long int));
    unsigned long long int *w0 = (unsigned long long int*)calloc(n_row2*m, sizeof(unsigned long long int));
    double *xw = (double*)calloc(10, sizeof(double));
    fread(a1, sizeof(unsigned long long int), n_row1,fp_a11);
    fread(a0, sizeof(unsigned long long int), n_row1,fp_a01);
    fread(w1, sizeof(unsigned long long int), n_row2,fp_a12);
    fread(w0, sizeof(unsigned long long int), n_row2,fp_a02);
//    printf("--------\n a1 %lld \n--------\n",a1[1]);
    
    //remove missing genotype
    for(i=0; i<n_row; i++){
        for(j=0; j<64; j++){
            temp1 = a0[i] & 1;
            temp2 = a1[i] & 1;
            temp3 = w0[i] & 1;
            temp4 = w1[i] & 1;
            temp5 = temp1+temp2;
            temp6 = temp3+temp4;
            if(temp5<2 && temp6<2){
                a0[n_row+i] <<= 1;
                a1[n_row+i] <<= 1;
                w0[n_row+i] <<= 1;
                w1[n_row+i] <<= 1;
                a0[n_row+i] = a0[n_row+i] | temp1;
                a1[n_row+i] = a1[n_row+i] | temp2;
                w0[n_row+i] = w0[n_row+i] | temp3;
                w1[n_row+i] = w1[n_row+i] | temp4;
                n_marker++;
            }
            a0[i] >>= 1;
            a1[i] >>= 1;
            w0[i] >>= 1;
            w1[i] >>= 1;
        }
    }
    get_n1234(&a1[n_row], &a0[n_row], &w1[n_row], &w0[n_row], n_row, xw);
    kin_robust(xw);
    printf("marker = %d\n",n_marker);
    fprintf(fp_r, "KING\t%.10lf\n", xw[5]);
    free(a1);
    free(a0);
    free(w1);
    free(w0);
    free(xw);
    fclose(fp_a11);
    fclose(fp_a01);
    fclose(fp_a12);
    fclose(fp_a02);
    fclose(fp_r);
    return 0;
}

int KING_single(char *input1, char *input2, char *output){
    FILE *fp_a1=NULL, *fp_a0=NULL, *fp_r=NULL;
    int m=2, i=0, j=0, n_row=0, n_marker=0;
    unsigned long long int temp1=0, temp2=0, temp3=0, temp4=0, temp5=0, temp6=0;
    char fileName_temp[1000] = "";
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".pre");
    if((fp_a1=fopen(fileName_temp,"rb"))==NULL){
        printf("The file .pre can't be opened!\n");
        return -1;
    }
    strcpy (fileName_temp,getFileName(input1));
    strcat (fileName_temp,getFileName(input2));
    strcat (fileName_temp,".pos");
    if((fp_a0=fopen(fileName_temp,"rb"))==NULL){
        printf("The file .pos can't be opened!\n");
        return -1;
    }
    strcpy (fileName_temp,output);
    strcat (fileName_temp,".king.txt");
    if((fp_r=fopen(fileName_temp,"a"))==NULL){
        printf("The .king.txt file isn't exist!\n");
        return -1;
    }
    fread(&n_row, sizeof(int), 1,fp_a0);
    fread(&n_marker, sizeof(int), 1,fp_a0);
    unsigned long long int *a1 = (unsigned long long int*)calloc(n_row*m, sizeof(unsigned long long int));
    unsigned long long int *a0 = (unsigned long long int*)calloc(n_row*m, sizeof(unsigned long long int));
    unsigned long long int *w1 = (unsigned long long int*)calloc(n_row*m, sizeof(unsigned long long int));
    unsigned long long int *w0 = (unsigned long long int*)calloc(n_row*m, sizeof(unsigned long long int));
    double *xw = (double*)calloc(10, sizeof(double));
    fread(a1, sizeof(unsigned long long int), n_row,fp_a1);
    fread(a0, sizeof(unsigned long long int), n_row,fp_a0);
    fread(w1, sizeof(unsigned long long int), n_row,fp_a1);
    fread(w0, sizeof(unsigned long long int), n_row,fp_a0);
//    printf("--------\n a1 %lld \n--------\n",a1[1]);
    //remove missing genotype
    n_marker=0;
    for(i=0; i<n_row; i++){
            for(j=0; j<64; j++){
                temp1 = a0[i] & 1;
                temp2 = a1[i] & 1;
                temp3 = w0[i] & 1;
                temp4 = w1[i] & 1;
                temp5 = temp1+temp2;
                temp6 = temp3+temp4;
                if(temp5<2 && temp6<2){
                    a0[n_row+i] <<= 1;
                    a1[n_row+i] <<= 1;
                    w0[n_row+i] <<= 1;
                    w1[n_row+i] <<= 1;
                    a0[n_row+i] = a0[n_row+i] | temp1;
                    a1[n_row+i] = a1[n_row+i] | temp2;
                    w0[n_row+i] = w0[n_row+i] | temp3;
                    w1[n_row+i] = w1[n_row+i] | temp4;
                    n_marker++;
                }
                a0[i] >>= 1;
                a1[i] >>= 1;
                w0[i] >>= 1;
                w1[i] >>= 1;
            }
        }
        get_n1234(&a1[n_row], &a0[n_row], &w1[n_row], &w0[n_row], n_row, xw);
        kin_robust(xw);
//        printf("king = %.4lf\n",xw[5]);
        fprintf(fp_r, "%d\t%.4lf\n", n_marker, xw[5]);
        free(a1);
        free(a0);
        free(w1);
        free(w0);
        free(xw);
        fclose(fp_a1);
        fclose(fp_a0);
        fclose(fp_r);
        return 0;
}

int single_run(char *fileName_input1, char *fileName_input2, char *path_input1, char *path_input2, const int type, char *fileName_output){
    char fileName_temp[10000];
//    char *temp_folder = malloc(1000);
    
//    char *input2[100] = {"bcftools", "merge", "--force-samples", "-g", "-", "/home/meng/KING/test_king/Test1_g1.b1.i1.bcf", "/home/meng/KING/test_king/Test1_g1.b1.i1.bcf", "-Oz", "-o", "/home/meng/KING/data_temp/Test1_g1.b1.i1.bcfTest1_g1.b1.i1.bcf.vcf.gz", (char*)NULL};
    
    FILE *fp_r=NULL;
    time_t start1, finish1;
    double duration1;
    start1 = time(0);
//    sprintf(temp_folder,"%s%d/",path_input2,(int)start1);//for windows change "/" to "\\"
//    strcpy (fileName_temp,"mkdir ");
//    strcat (fileName_temp,temp_folder);
//    printf("%s\n",fileName_temp);
//    system(fileName_temp);
    strcpy (fileName_temp,fileName_output);
    strcat (fileName_temp,"king_chi.log");
    if((fp_r=fopen(fileName_temp,"w"))==NULL){
        printf("Couldn't creat the log file!\n");
        return -1;
    }
    
    if(type==1){
        strcpy (fileName_temp,path_input1);
        strcat (fileName_temp,"bcftools merge --force-samples -g - -R ");
        
        strcat (fileName_temp,path_input1);
        strcat (fileName_temp,"hg38.sites2include.bed ");
        strcat (fileName_temp,fileName_input1);
        strcat (fileName_temp," ");
        strcat (fileName_temp,fileName_input2);
        strcat (fileName_temp," -Ob -o temp.bcf");
    //    strcat (fileName_temp," -Ou | bcftools filter -g 5 -e \"FORMAT/DP<1\" |  bcftools view -v snps -f PASS -Ov | python3 "); // for windows
        
        system(fileName_temp);
        fprintf(fp_r, "%s\n", fileName_temp);
    }
    
    strcpy (fileName_temp,path_input1);
    strcat (fileName_temp,"plink2 --bcf temp.bcf --make-bed  --snps-only --autosome --max-alleles 2 --set-all-var-ids @:#:\\$r:\\$a --out temp > temp.plink2_std");
    system(fileName_temp);
    fprintf(fp_r, "%s\n", fileName_temp);

    strcpy (fileName_temp,path_input1);
    strcat (fileName_temp,"plink2 --bfile temp --export A --out temp >> temp.plink2_std");
    system(fileName_temp);
    fprintf(fp_r, "%s\n", fileName_temp);


    if (compress_bit(fileName_input1, fileName_input2, fileName_output)!=0) {
        printf("\nThe compressed binary file %s and %s can't be created!\n", fileName_input1, fileName_input2);
        return -1;
    }

    if (KING_single(fileName_input1, fileName_input2, temp_folder, fileName_output)!=0) {
        printf("\nThe King coefficient calculation of file %s and %s can't be finished!\n", fileName_input1, fileName_input2);
        return -1;
    }
    finish1 = time(0);
    duration1 = difftime(finish1, start1);
//    printf( "Time consumption = %.0f seconds\n", duration1);
    fprintf(fp_r, "Time consumption = %.0f seconds\n", duration1);
    fclose(fp_r);
    return 0;
}

int single_run2(char *fileName_input1, char *fileName_input2, char *path_input1, char *path_input2, char *fileName_output){
    char fileName_temp[1000];
    char *temp_folder = malloc(1000);
//    char *input2[100] = {"bcftools", "merge", "--force-samples", "-g", "-", "/home/meng/KING/test_king/Test1_g1.b1.i1.bcf", "/home/meng/KING/test_king/Test1_g1.b1.i1.bcf", "-Oz", "-o", "/home/meng/KING/data_temp/Test1_g1.b1.i1.bcfTest1_g1.b1.i1.bcf.vcf.gz", (char*)NULL};
    FILE *fp_r=NULL;
    time_t start1, finish1;
    double duration1;
    start1 = time(0);
    sprintf(temp_folder,"%s%d/",path_input2,(int)start1);//for windows change "/" to "\\"
    strcpy (fileName_temp,"mkdir ");
    strcat (fileName_temp,temp_folder);
//    printf("%s\n",fileName_temp);
    system(fileName_temp);
    strcpy (fileName_temp,temp_folder);
    strcat (fileName_temp,"king_chi.log");
    if((fp_r=fopen(fileName_temp,"w"))==NULL){
        printf("Couldn't creat the log file!\n");
        return -1;
    }
    
    strcpy (fileName_temp,path_input1);
    strcat (fileName_temp,"plink2 --vcf ");
    strcat (fileName_temp,fileName_input1);
    strcat (fileName_temp," --make-bed  --snps-only --autosome --max-alleles 2 --set-all-var-ids @:#:\\$r:\\$a --out ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp," > ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,".plink2_std");
    system(fileName_temp);
    fprintf(fp_r, "%s\n", fileName_temp);
        
    strcpy (fileName_temp,path_input1);
    strcat (fileName_temp,"plink2 --vcf ");
    strcat (fileName_temp,fileName_input2);
    strcat (fileName_temp," --make-bed  --snps-only --autosome --max-alleles 2 --set-all-var-ids @:#:\\$r:\\$a --out ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp," > ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".plink2_std");
    system(fileName_temp);
    fprintf(fp_r, "%s\n", fileName_temp);
    
    if (file_check(fileName_input1, fileName_input2, temp_folder, fileName_output)<0) {
        printf("\nPlease check input file %s and %s\n", fileName_input1, fileName_input2);
        return -1;
    }
    
    strcpy (fileName_temp,path_input1);
    strcat (fileName_temp,"plink2 --bfile ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp," --extract ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_output));
    strcat (fileName_temp,".snp_list.txt --export A --out ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp," >> ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input1));
    strcat (fileName_temp,".plink2_std");
    system(fileName_temp);
    fprintf(fp_r, "%s\n", fileName_temp);
    
    strcpy (fileName_temp,path_input1);
    strcat (fileName_temp,"plink2 --bfile ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp," --extract ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_output));
    strcat (fileName_temp,".snp_list.txt --export A --out ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp," >> ");
    strcat (fileName_temp,temp_folder);
    strcat (fileName_temp,getFileName(fileName_input2));
    strcat (fileName_temp,".plink2_std");
    system(fileName_temp);
    fprintf(fp_r, "%s\n", fileName_temp);
    
    
    if (compress_bit2(fileName_input1, temp_folder, fileName_output)!=0) {
        printf("\nThe compressed binary file %s can't be created!\n", fileName_input1);
        return -1;
    }
    
    if (compress_bit2(fileName_input2, temp_folder, fileName_output)!=0) {
        printf("\nThe compressed binary file %s can't be created!\n", fileName_input2);
        return -1;
    }

    if (KING(fileName_input1, fileName_input2, temp_folder, fileName_output)!=0) {
        printf("\nThe King coefficient calculation of file %s and %s can't be finished!\n", fileName_input1, fileName_input2);
        return -1;
    }
    finish1 = time(0);
    duration1 = difftime(finish1, start1);
//    printf( "Time consumption = %.0f seconds\n", duration1);
    fprintf(fp_r, "Time consumption = %.0f seconds\n", duration1);
    fclose(fp_r);
    
    remove_temp(fileName_input1, fileName_input2, temp_folder);
    strcpy (fileName_temp,temp_folder);
    remove(fileName_temp);
    strcpy (fileName_temp,"rm -r ");
//    strcpy (fileName_temp,"rmdir /s /q "); //for windows
    strcat (fileName_temp,temp_folder);
    system(fileName_temp);
    free(temp_folder);
    return 0;
}

int multiple_run(char *fileName_input1, char *fileName_input2, char *path_input1, char *path_input2, const int type, char *fileName_output){
    int n=0;
//    char *input1 = "";
    char *input2 = "";
    FILE *fp_bim2=NULL;
    char fileName_temp[100] = "";
//    strcpy (fileName_temp,fileName_input1);
//    if((fp_bim1=fopen(fileName_temp,"r"))==NULL){
//        printf("The file list1 doesn't exist!\n");
//        return -1;
//    }
    strcpy (fileName_temp,fileName_input2);
    if((fp_bim2=fopen(fileName_temp,"r"))==NULL){
        printf("The file list2 doesn't exist!\n");
        return -1;
    }
    
    struct gm{
        char name[1000];
    };
//    struct gm a[2]={0};
    struct gm b[2]={0};
//    if(fscanf(fp_bim1,"%s", a[0].name)>0){
//        input1 = a[0].name;
//    }else{
//        printf("The file list1 is empty!\n");
//        return -1;
//    }
    while(fscanf(fp_bim2,"%s", b[0].name)>0){
        input2 = b[0].name;
        n++;
        single_run(fileName_input1, input2, path_input1, path_input2, type, fileName_output);
//        printf("The file list2 is %s %s\n",input1, input2);
    }
    if(n==0){
        printf("The file list2 is empty!\n");
        return -1;
    }
    fclose(fp_bim2);
    return 0;
}



int main(int argc, char * argv[]) {
    char *fileName_input1 = "/eva/projects/gengen/ShareTest/data/7035.7046.bcf.vcf";//target file name
    char *fileName_input2 = "/eva/projects/gengen/ShareTest/data/sample1.bcf";//data base file name
	char *path_input1 = "./";//pathway of app
    char *path_input2 = "./";//pathway of temp files
    char *fileName_output = "myKING";
    char fileName_temp[1000];
    int i=0, type=2, run=3;
    FILE *fp_r=NULL;
    char *s1, *s2, *s3, *s4, *s5, *s6, *s7;
    s1 = "--file";
    s2 = "--data";
    s3 = "--out";
    s4 = "--format";// type=1 means bcf(gvcf); type=2 means vcf;
    s5 = "--module";//run=2 means loop the file list one by one; run=1 means conduct two input file only; run=3 means conduct two input vcf file only;
    s6 = "--path_app";
    s7 = "--path_temp";
    
    for (i = 0; i < argc; i++) {
        if (strncmp(s1, argv[i], 6)==0) {
            fileName_input1 = argv[i+1];
//            printf(" --file %s\n", fileName_input1);
            fileName_output= argv[i+1];
        }
        if (strncmp(s2, argv[i], 6)==0) {
            fileName_input2 = argv[i+1];
//            printf(" --data %s\n", fileName_input2);
            fileName_output= argv[i+1];
        }
        if (strncmp(s3, argv[i], 5)==0) {
            fileName_output = argv[i+1];
//            printf(" --out %s\n", fileName_output);
        }
        if (strncmp(s4, argv[i], 8)==0) {
            type = atoi(argv[i+1]);
//            printf(" --format %d\n", type);
        }
        if (strncmp(s5, argv[i], 8)==0) {
            run = atoi(argv[i+1]);
//            printf(" --module %d\n", run);
        }
        if (strncmp(s6, argv[i], 10)==0) {
            path_input1 = argv[i+1];
//            printf(" --path_app %s\n", path_input1);
        }
        if (strncmp(s7, argv[i], 11)==0) {
            path_input2 = argv[i+1];
//            printf(" --path_temp %s\n", path_input2);
        }
    }
    
    
    strcpy (fileName_temp,fileName_output);
    strcat (fileName_temp,".king.txt");
    if((fp_r=fopen(fileName_temp,"w"))==NULL){
        printf("The .king.txt file isn't exist!\n");
        return -1;
    }
    fprintf(fp_r, "#Sample_ID\tReference_ID\tNumber_of_markers\tKinship_coefficient\n");
    fclose(fp_r);
    
    if (run == 1){
        if (single_run(fileName_input1, fileName_input2, path_input1, path_input2, type, fileName_output)!=0) {
            printf("\nThe King coefficient calculation can't be finished!\n");
            return -1;
        }
    }
    
    
    }
    
    
    return 0;
}







