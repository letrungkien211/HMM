hmm
===

Hidden markov model.

From "Introduction to machine learning" by Ethem Alpaydin
Theories:

***********************************
Arguments:
1. N: Number of states in the model
   S = {S(1), S(2),...,S(n)}
2. M: Number of distinct observation symbols in the alphabet
   V = {v(1), v(2), ..., v(m)}
3. State transition probabilities
   A = {a(i,j)} where a(i,j) = P(q(t+1) = Sj | q(t) = S(i))
4. Observation probabilities
   B = [b(j,m)] where b(j,m) = P(O(t) = v(m)| q(t) = S(j)
5. Initial state probabilities
   pi = [pi(i)] where pi(i) = P(q(1) = Si)

***********************************
Three basic problems of HMMs
Given a number of sequences of obsevervations:
1. Given a model Lambda, we would like to evaluate the probability of any give observation sequence, O , namely, P(O|Lambda)
2. given a model Lambda and an observation O, we would like to find out the state sequence Q, which has the highest probability of generating O; namely, we want to find Q* that maximizes P(Q|O, Lambda)
3. Given a training set of observation sequences, X = {O(k)}(k=1,..,K}, we would like to learn the model that maximizes the probability of generating X; namely, we want to find Lambda* that maximizes P(X|Lambda)


***********************************
