#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#include "../costMatrix_2d.hpp"
#include "../dynamicCharacterOperations.h"

#define __STDC_FORMAT_MACROS

// #include "seqAlignForHaskell.h"

int main() {
    const size_t alphabetSize = 5;
    const size_t tcmLen       = alphabetSize * alphabetSize;

    unsigned int tcm[tcmLen];
    for (size_t i{0}; i < alphabetSize; i++) {
        for (size_t j{0}; j < alphabetSize; j++) {
            if (i == j) {
                tcm[i * alphabetSize + j] = 0;
            }
            else {
                tcm[i * alphabetSize + j] = 1;
            }
        }
    }



    /*
    const size_t SEQ_A_LEN = 15;
    packedChar seqA_main[SEQ_A_LEN];

    for (size_t i = 0; i < SEQ_A_LEN; i++) {
        seqA_main[i] = 15 - i;
    }

    const size_t SEQ_B_LEN = 10;
    packedChar seqB_main[SEQ_B_LEN];
    for (size_t i = 0; i < SEQ_B_LEN; i++) {
        seqB_main[i] = i;
    }
    */


    CostMatrix_2d myMatrix = CostMatrix_2d(alphabetSize, tcm);


    auto firstKey  = makeDCElement( alphabetSize, 1 );
    auto secondKey = makeDCElement( alphabetSize, 1 );
    auto retMedian = makeDCElement( alphabetSize, 1 );
    auto foundCost{0};

    // just a test: alphabet size == 4, so don't need packedChar*
    auto median = CANONICAL_ZERO;

    // First, test constructor, i.e. that unambiguous characters have been inserted.
    printf("\n\n\n******* Testing constructor: insertion of unambiguous characters. ******\n");
    for (size_t key1 = 0; key1 < alphabetSize; ++key1) { // for every possible value of key1, key2
        SetBit(firstKey->element, key1);
        SetBit(&median, key1);    // computed median just for testing.
        // printPackedChar(&median, 1, alphabetSize);

        for (size_t key2 = 0; key2 < alphabetSize; ++key2) { // no longer assumes 0 diagonal
            SetBit(secondKey->element, key2);
            //auto cost = tcm[key1 * alphabetSize + key2];
            SetBit(&median, key2);

            foundCost = myMatrix.getSetCostMedian(firstKey, secondKey, retMedian);
            fflush(stdout);
            printf("key 1 set: %zu\n", key1);
            printf("key 2 set: %zu\n", key2);
            printf("found median:\n");
            printPackedChar(retMedian->element, 1, alphabetSize);
            printf("found cost:    %d\n", foundCost);

            if(key2 != key1) ClearBit(&median, key2); // the key1 bit needs to persist on the median
            ClearBit(secondKey->element, key2);
        } // key2
        ClearBit(firstKey->element, key1);
        ClearBit(&median, key1);
    }
    printf("Passed!\n\n\n");

    printf("\n\n\n******* Testing ambiguous characters: get/set of ambiguous characters. ******\n");
    size_t numSetInKey;
    for(size_t i = 0; i < 25; ++i) {
        printf("\n\niteration %2zu\n", i + 1);
        numSetInKey = rand() % alphabetSize + 1;

        for(size_t setIdx = 0; setIdx < numSetInKey; ++setIdx) {
            SetBit(firstKey->element, rand() % alphabetSize);
        }

        numSetInKey = rand() % alphabetSize + 1;
        for(size_t setIdx = 0; setIdx < numSetInKey; ++setIdx) {
            SetBit(secondKey->element, rand() % alphabetSize);
        }

        printf("Input:\n");
        printf("  Key #1: "); printPackedChar(firstKey->element,  1, alphabetSize); printf("\n");
        printf("  Key #2: "); printPackedChar(secondKey->element, 1, alphabetSize); printf("\n");
        foundCost = myMatrix.getSetCostMedian(firstKey, secondKey, retMedian);
        printf("Output:\n");
        printf("  Median: "); printPackedChar(retMedian->element, 1, alphabetSize); printf("\n");
        printf("  Cost: %i\n", foundCost);

        ClearAll( firstKey->element, dynCharSize(alphabetSize, 1) );
        ClearAll(secondKey->element, dynCharSize(alphabetSize, 1) );
    }

    // Free everything we have allocated as to not mess with valgrind's leak diagnostics.
    freeDCElem(firstKey);
    freeDCElem(secondKey);
    freeDCElem(retMedian);
    free(firstKey);
    free(secondKey);
    free(retMedian);


    /****** An abandoned packedCharOr test that doesn't scale with alphabetSize > 64 ******/

    /**
    auto first   = new packedChar{1,0},
         second  = new packedChar{4,0},
	 third   = new packedChar{16,0},
         result  = packedCharOr( first, second, alphabetSize, 1 ),
         result2 = packedCharOr( result, third, alphabetSize, 1 );


    printf("%" PRIu64 "\n", *result);
    printf("%" PRIu64 "\n\nDone.", *result2);

    delete first;
    delete second;
    delete third;
    free(result);
    free(result2);
    **/


    /****** This next to test Yu Xiang's code, once you can. ******/

    // int success = aligner(seqA_main, SEQ_A_LEN, seqB_main, SEQ_B_LEN, alphabetSize, getCostMatrix(myMatrix), &retMedChar);

    // if (success == 0) {
    //     printf("\nSuccess!\n\n");
    //     // printf("The aligned sequences are: \n%p\n%p\n", retAlign->seq1, retAlign->seq2);
    //     // printf("The cost of the alignment is: %d\n", retAlign->weight);
    //     // for(int i = 0; i < length; ++i) {
    //     //     printf("%d\n",(int)retAlign->seq1[i]);
    //     // }
    // } else {
    //     printf("Fail!\n");
    // }

    // free(seqA_main);
    // free(seqB_main);

}
