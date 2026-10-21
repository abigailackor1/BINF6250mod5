# Introduction

---

## Smith–Waterman: Local Sequence Alignment

The Smith–Waterman algorithm is a dynamic programming algorithm used to identify the optimal **local alignment** between two DNA, RNA, or protein sequences. Unlike global alignment, which compares sequences across their entire lengths, local alignment identifies the highest-scoring regions within the sequences.

Dynamic programming solves a complex problem by dividing it into smaller, overlapping subproblems and storing their solutions. In the Smith–Waterman algorithm, these solutions are stored in a scoring matrix.

## Algorithm

1. Define a scoring scheme, such as:
   - Match: `+1`
   - Mismatch: `-1`
   - Gap: `-2`

2. Place the first sequence horizontally across the top of the scoring matrix.

3. Place the second sequence vertically along the left side of the matrix.

4. Create a scoring matrix with one additional row and column for initialization. If the sequences have lengths `m` and `n`, the matrix dimensions are `(m + 1) × (n + 1)`.

5. Initialize the first row and first column of the matrix to `0`.

6. Calculate the score for each remaining cell using four possible values:

   - **Diagonal:** previous diagonal value + match or mismatch score
   - **Up:** value above the current cell + gap penalty
   - **Left:** value to the left of the current cell + gap penalty
   - **Zero:** `0`

7. Assign the maximum of these four values to the current cell:

   `H(i, j) = max(0, diagonal score, up score, left score)`

   Including zero prevents negative values from appearing in the matrix and allows a new local alignment to begin at any position.

8. Repeat the calculation until the entire scoring matrix is filled.

9. Begin the traceback at the highest-scoring cell in the matrix.

10. Follow the path that produced each score:

    - A diagonal move aligns two sequence characters.
    - An upward move introduces a gap in the horizontal sequence.
    - A leftward move introduces a gap in the vertical sequence.

11. End the traceback when a cell containing `0` is reached. The resulting path represents an optimal local alignment under the selected scoring scheme.

## Notebook Overview

This notebook implements the Smith–Waterman algorithm and evaluates it using simulated sequence data.

Three test cases are included:

1. **Reference example:** Reproduces the sequence-alignment example provided in the lecture slides.
2. **Repeated regions:** Evaluates whether the algorithm can identify locally similar regions within sequences containing repeated patterns.
3. **Negative-score handling:** Confirms that negative candidate scores are replaced with zero, as required by the Smith–Waterman algorithm.

The implementation begins with helper functions that calculate the diagonal, up, and left scores. Separating these calculations into helper functions improves readability and keeps each function focused on a single task.

The remaining functions construct and fill the scoring matrix, perform the traceback, and return the optimal local alignment produced by the Smith–Waterman algorithm.

# Pseudocode

---

Put pseudocode in this box:

```
Helper Functions: CALCULATE the scores for above, left and diagonal 

def get_node_above()
    GET value from the cell above
    IF greater than 0
        ADD either match score, mismatch penalty or gap penalty
    ELSE return 0
  
def get_node_left()
    GET value from the cell to the left
    IF greater than 0
        ADD either match score, mismatch penalty or gap penalty
    ELSE return 0

def get_node_diag()
    GET value from the cell diagonal
    IF greater than 0
        ADD either match score, mismatch penalty or gap penalty
    ELSE return 0
```
```
Main functions: Perform Smith-Waterman algorithm for optimal aligned sequences

def cal_score()
    left ←  CALL get_node_left(value to the left of matrix[i][j])
    up ← CALL get_node_above(value above matrix[i][j])
    diagonal ← CALL get_node_diag(value diagonally above-left of matrix[i][j])

    IF seq1[i - 1] equals seq2[j - 1]:
        diagonal_score ← diagonal + match

    ELSE:
        diagonal_score ← diagonal + mismatch

    up_score ← up + gap penalty
    left_score ← left + gap penalty

    score ← MAX(0, diagonal_score, up_score, left_score)
    matrix[i][j] ← score

    IF score equals 0:
        trace ← END
    ELSE IF score equals diagonal_score:
        trace ← DIAGONAL
    ELSE IF score equals up_score:
        trace ← UP
    ELSE:
        trace ← LEFT

    RETURN score, trace

def traceback()
    aligned_seq1 ← empty string
    aligned_seq2 ← empty string

    (i, j) ← maximum_position

    WHILE i > 0 AND j > 0:
        current_move ← traceback_matrix[i][j]

        IF current_move = DIAGONAL:
            APPEND seq1[i - 1] to aligned_seq1
            APPEND seq2[j - 1] to aligned_seq2
            MOVE diagonally up-left

        ELSE IF current_move = UP:
            APPEND seq1[i - 1] to aligned_seq1
            APPEND a gap to aligned_seq2
            MOVE up

        ELSE IF current_move = LEFT:
            APPEND a gap to aligned_seq1
            APPEND seq2[j - 1] to aligned_seq2
            MOVE left

        ELSE:
            STOP traceback

    REVERSE aligned_seq1
    REVERSE aligned_seq2

    RETURN aligned_seq1, aligned_seq2

def smith_waterman()
    CREATE score_matrix and traceback_matrix with dimensions:
        (length of seq1 + 1) × (length of seq2 + 1)

    INITIALIZE both matrices with zeros

    SET match = 1
    SET mismatch = -1
    SET gap = -1
    SET max_score = 0

    FOR each position i in seq1:
        FOR each position j in seq2:
            (score, trace) ← CalculateScore(
                score_matrix, seq1, seq2,
                i, j, match, mismatch, gap
            )

            STORE score in score_matrix[i][j]
            STORE trace in traceback_matrix[i][j]

            Check that score is not negative

            IF score is greater than max_score:
                max_score ← score

    FIND all positions in score_matrix equal to max_score

    FOR each maximum-score position:
        PERFORM Traceback from that position
        STORE the resulting alignment

    SELECT the first stored alignment

    RETURN aligned_seq1, aligned_seq2, and score_matrix

```

# Biological Considerations and Context 

--- 

Several approaches are available for sequence alignment, and each is suited to a different biological question:

- **Global alignment:** Needleman–Wunsch
- **Local alignment:** Smith–Waterman
- **Multiple sequence alignment:** MUSCLE
- **Heuristic sequence search:** BLAST

The Smith–Waterman algorithm performs pairwise local alignment. It identifies the highest-scoring region shared by two sequences rather than aligning them across their entire lengths. Although it guarantees an optimal local alignment for the selected scoring system, several computational and biological limitations should be considered.

## Limitations and Testing Considerations

### 1. Scalability

**Limitation:**  
The standard Smith–Waterman algorithm compares two sequences at a time. It is therefore not a multiple sequence alignment method. Its computational requirements also increase as the sequences become longer, making direct comparisons against large sequence databases expensive.

This limitation does not impose a fixed maximum sequence length. Instead, the practical limit depends on factors such as the available memory, processing time, implementation, and computing hardware.

**Possible tests:**

- Measure execution time as sequence length increases.
- Measure memory usage for different matrix dimensions.
- Test sequences of increasing length to determine the practical limits of this implementation.
- Compare the results and runtime with a heuristic search tool such as BLAST.

### 2. Sensitivity to the Scoring Scheme

**Limitation:**  
The selected match score, mismatch penalty, and gap penalty directly affect the scoring matrix and traceback path. A poorly chosen scoring scheme may produce an alignment that is mathematically optimal but biologically uninformative.

For protein sequences, substitution matrices such as BLOSUM or PAM are generally more informative than assigning the same score to every amino acid match or mismatch. Gap penalties should also reflect the biological assumption that insertions and deletions may occur.

**Possible tests:**

- Run the same sequence pair with several scoring schemes.
- Compare how different mismatch and gap penalties affect the aligned region.
- Determine whether changing the scoring scheme alters the maximum score or traceback path.

### 3. Time and Space Complexity

**Limitation:**  
For sequences of lengths `m` and `n`, the standard Smith–Waterman implementation has:

- **Time complexity:** `O(mn)`
- **Space complexity:** `O(mn)`

Every cell in an `(m + 1) × (n + 1)` matrix must be calculated and stored when the complete matrix is retained for traceback. As a result, execution time and memory use increase approximately with the product of the sequence lengths.

**Possible tests:**

- Record runtime and memory use for sequences of increasing length.
- Compare observed performance with the expected `O(mn)` growth.
- Test memory-efficient implementations when the full scoring matrix is not required.
- Use heuristic tools such as BLAST when searching large databases, while recognizing that a heuristic method may not guarantee the same optimal alignment as Smith–Waterman.

### 4. Sequence Order and Orientation

**Limitation:**  
A standard Smith–Waterman alignment assumes that the aligned residues occur in the same relative order and orientation within each aligned region. A single alignment may not adequately represent complex genomic events such as inversions, translocations, or large duplications.

Separate local alignments may detect related regions, but additional analysis is often required to determine their genomic order, orientation, and broader structural relationship.

**Possible tests:**

- Align each DNA sequence against both the forward sequence and its reverse complement.
- Introduce simulated inversions, translocations, and duplications to evaluate how the implementation responds.
- Examine multiple high-scoring local alignments instead of reporting only one optimal alignment.
- Compare the results with tools designed for genome alignment or structural-variation analysis.

## Biological Interpretation

A high Smith–Waterman score represents the best local alignment under the selected scoring scheme. However, it does not by itself establish homology, shared function, or biological significance. Additional evidence, such as statistical significance, sequence annotations, conserved domains, genomic context, and evolutionary relationships, should be considered.

In practice, bioinformatics analysis follows a connected and often iterative pipeline. Researchers do not necessarily complete one step, accept its output without review, and permanently move on. Instead, the results of later stages may reveal problems that require revisiting data preparation, parameter selection, alignment, or another earlier stage.

Quality-control checks should therefore be performed before and after each major stage. Before alignment, the sequence data should be checked for issues such as incorrect formatting, ambiguous characters, contamination, or incomplete sequences. After alignment, the scores, aligned regions, gaps, and biological plausibility of the results should be evaluated before the data are used in downstream analyses. This process limits the propagation of errors and supports reliable, reproducible biological conclusions.

## In Comparison to EMBOSS Water 

The EMBOSS Water tool utilizes the Smith-Waterman algorithm with modified parameters for speed enhancements. Out of curiosity, we used the test sequence given from the lecture slides to compare our code to an actual standardized tool. Here are the parameters we used on the tool interface: 
1. Sequence Type: DNA
2. First Sequence:

>test1
TACTTAG

3. Second Sequence:

>test2
CACATTAA

5. Output format: Pair (to perform pairwise alignment)
6. Title: Test

Here is the outcome: 

<img width="628" height="672" alt="Screenshot 2026-10-21 at 2 18 47 AM" src="https://github.com/user-attachments/assets/d00ef578-ebae-4fb5-9238-d0fffeef2321" />


### Mini-Reflection on Our Code and EMBOSS Water

It is really interesting to view the differences in alignment. As we had mentioned before, the scoring scheme does have an impact on alignment, as our results produce a different output than EMBOSS Water. If we had more time, we would have taken additional steps to learn more about the tool and its outputs along with a comparison to other existing alignment tools. We think that this highlights the issue of standardization and reproducibility in Bioinformatics; with many new tools and researchers creating their own pipelines, it can be difficult to traceback results due to the various options of tools and parameters. 

# Successes
**Matrix traversal:** We were able to understand the premise of the project, and successfully store and access scores is terms of `match`, `mistmatch`, and `gap` points in `score_matrix`. Additionally, we were able to store directions in our `traceback_matrix` (`left`, `up`, and `diag`) and use those to trace back from `max_pos` where we stored the maximum score of the alignment.

**Scoring strategy:** We were able to successfully implement a dynamic scoring system that allowed for customizable match, mistmatch, and gap penalties.

**Biological understanding:** We were able to grasp an understanding of how the Smith-Waterman algorithm is applied to optimize the local alignment of two sequences. We learned that there will be limitations to the overall algorithm and our code due to variability. 

**Comparison to an Actual Tool:** We were able to set aside time to try and compare the outputs of our code to an existing tool. This validated our limitation that the scoring scheme truly does matter and can affect alignment. One has to consider how we would score matches, mismatches and gaps and determine the weight of each one. It is not a simple "choose whichever numbers we would like to pick". 

# Struggles
**Scoring parameters:** Initially, we had some confusion when it came to the scoring. We misinterpreted the instruction and comments from the class, thinking that we should not include pentalies with negative numbers. It took playing around with the scoring metrics and conducting a bit of further research, to realize that it was in fact the scores themselves that were stored in the matrix that should not be negative. We addressed this by adding a one line check to make sure that if the score in the current matrix was negative, we should replace it with `0`.

**Reversal of strands:** We also had to print out our aligned sequences to realize that they were returning in reverse order. This was also an easy fix, by using `aligned_seq1 = aligned_seq1[::-1]` and `aligned_seq2 = aligned_seq2[::-1]`

**Numpy Broadcasting**: Due to a time constraint, we were unable to implement the broadcasting function in Numpy. From our understanding this function would "arrays with different shapes during arithmetic operations" as said from the Numpy manual. We learned that this function can be helpful in algorithm implementation, but also has its constraints as well (e.g. inefficient use of memory that can slow computation). 

# Personal Reflections
## Group Leader: Maggie Wenger
This project was a lot of fun for me! I think it was great to have a project that was a little bit easier to understand from the jump, and allowed me to focus more heavily on the biological context and usage for this algorithm. I think it was also great practice for me with matrices: I have had exposure to them and coded with them before, but not for a while, so it was nice to dive back into that previous knowledge. 

Something I would have like to implement if we had more time would be the use of numpy broadcasting. I remember we spoke briefly about it during class, and upon research seems like a very efficient way to edit and apply calculations to matrices. I also spent a lot of time playing around with different scoring options, as it seems biologically that a gap should be penalized more than a mismatch. I couldn't come up with a way to make this work functionally, but would have liked some more time to play around with these parameters.

## Other member
Other members' reflections on the project
Abby - I really enjoyed this project, i feel like it allowed me to slow down after the past few weeks and really focus on the more meticulous parts of my problem solving process. The methodology made sense and having to dig deeper into the repercussions of biological aspects was also very interesting for this project. 

## Other member: Kailey Solorzano

The project was a refresher in comparison to previous projects. A learning curve was understanding the biological context over a computed output. I feel that this week, I was able to dive deep into the Smith-Waterman algorithm and take into consideration how this is used in genome alignment. The provided data is only a sample that was generated so I wonder how our code would compare to an actual targeted sequence. 

With this given project, I was able to take a moment to reflect on dynammic programming and using matrices on Python. My experience with matrices stems from using on R so it was interesting to utilize a different coding language for something similar. I also took a moment to reflect on dynamic programming as a whole- taking a large problem and creating smaller sub-problems that will make it easier to code, also that it is okay to not produce an output. This project altered my thinking on how to approach future projects as I learned to take a few days to digest the material and remember that bioinformatics is both biology and coding. 

# Generative AI Appendix
Anthropic, Sonnet 5.5. 2026. Claude [Large language model]. https://claude.ai

Prompt: Check the README for repetitive and confusing language. 

Output: As seen on the current README 

Justification: Often, I will word Introductions and other components of the README in a way that can be confusing to reviewers. The use of generative AI checks for my grammar issues and proofreads for repetitive language. 
