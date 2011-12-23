#include <iostream>
#include <vector>
#include <algorithm>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Utils.h>
#include <set>




using namespace std;


class Pair {// A pair will be defined by two genes and the score Delta and Gamma.
    public:
        double delta, gamma;
        int geneA, geneB;
        Pair (double, double,int, int);
};

Pair :: Pair (double d, double g, int a, int b){
    delta = d;
    gamma = g;
    geneA = a;
    geneB = b;
}

bool compare (Pair pair1, Pair pair2) {//Define a criteria to order two pairs based on theirs scores Delta and Gamma.
    if((pair1.delta > pair2.delta)|| (pair1.delta == pair2.delta && pair1.gamma > pair2.gamma)){
        return 1;
    }
    else return 0;
}



extern "C" {

SEXP cktspair(SEXP Rdat, SEXP Rgrp, SEXP Rnumberk, SEXP Rreplacena, SEXP Rlength)
{
    int i,j,k,u;
    int size;
    double score, rank;
    int nProtected = 0;
    int probsum1[2]; probsum1[0] = 0; probsum1[1] = 0;
    int probsum2[2]; probsum2[0] = 0; probsum2[1] = 0;
    int probrank[2]; probrank[0] = 0; probrank[1] = 0;

    double *grp, *dat, *repna;
    grp = REAL(Rgrp);
    dat = REAL(Rdat);
    repna = REAL(Rreplacena);

    int *length;
    length = INTEGER(Rlength);

    int *numberk;
    numberk = INTEGER(Rnumberk);

    int n, m;
    n = LENGTH(Rgrp);
    m = LENGTH(Rdat)/n;

    double n0, n1;
    n0 = 0;
    n1 = 0;

    vector<Pair> listPair;

    for(i = 0; i < m; i++){
        for(j = 0; j < i; j++){
            for(k = 0; k < n; k++){
		//Compute the scores for the pair (j,i)
                if(dat[i + k * m]!= *repna && dat[j + k * m]!= *repna){
                    if(grp[k] == 1){
                        probsum1[1] += (dat[i + k * m] < dat[j + k * m]);
                        probsum2[1] += (dat[i + k * m] > dat[j + k * m]);
                        probrank[1] += (dat[i + k *m] - dat[j + k * m]);
                        n1++;
                    }
                    if(grp[k] == 0){
                        probsum1[0] += (dat[i + k * m] < dat[j + k * m]);
                        probsum2[0] += (dat[i + k * m] > dat[j + k * m]);
                        probrank[0] += (dat[i + k *m] - dat[j + k * m]);
                        n0++;
                    }
                }
            }
	    if(n0 !=0 && n1 != 0){
	            score = max(fabs(probsum1[1]/n1 - probsum1[0]/n0),fabs(probsum2[1]/n1 - probsum2[0]/n0));
        	    rank = fabs(probrank[1]/n1 - probrank[0]/n0);
	    }	
            Pair pair (score, rank, j+1, i+1);
            size = listPair.size();

            if(size < length[0]){listPair.push_back(pair);}//Fill the list of the pair if it has not the desired length.

            if(size == length[0]){//Check if the pair (j,i) is better than the worst pair in the listPair.
		     if(compare(listPair[length[0]-1],pair)==0){//If it is better, it will replace the worst one.
       		 	 listPair.erase(listPair.begin() + length[0]-1);
       		 	 listPair.push_back(pair);
      			 sort(listPair.begin(), listPair.end(), compare);
    		     }
	   }

	    score = 0;
	    rank = 0;
            probsum1[0] = 0;
            probsum1[1] = 0;
            probsum2[0] = 0;
            probsum2[1] = 0;
            probrank[0] = 0;
            probrank[1] = 0;
            n0 = 0;
            n1 = 0;
            void R_CheckUserInterrupt(void);
        }
    }

	set <double> indices;
	set<double>::iterator itA, itB;

	SEXP rfinal;
	double *final;
	PROTECT(rfinal = allocVector(REALSXP,4 * numberk[0]));
	final = REAL(rfinal);
	nProtected++;
	size = 0;
	u=0;

	for(i=0; i <4*numberk[0]; i++){
		final[i]=0;
	}

	//Select only the numberk best pairs in listPair with the restriction that a gene can be used in only one pair.

        while((size < 2*numberk[0]) & (u < length[0])){
		itA = indices.find(listPair[u].geneA);
		itB = indices.find(listPair[u].geneB);
		if((itA == indices.end()) & (itB == indices.end())){
			size = indices.size();
			final[2*size] = listPair[u].delta;
			final[2*size+1] = listPair[u].gamma;
			final[2*size+2] = listPair[u].geneA;
			final[2*size+3] = listPair[u].geneB;
			indices.insert(listPair[u].geneA);
			indices.insert(listPair[u].geneB);
		}
		u++;
		size = indices.size();
        }


SEXP list;

PROTECT(list=allocVector(VECSXP,1));
++nProtected;
SET_VECTOR_ELT(list, 0, rfinal);
UNPROTECT(nProtected);
return(list);
}
}
