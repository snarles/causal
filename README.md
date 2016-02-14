# causal
Causal Inference stuff

# Abstract

Causal Inference and Causal Invariance

Part I: Introduction

by Charles Zheng and Qingyuan Zhao

In this two-part talk, we will give an overview to the field of causal inference, with the goal of discussing recent developments in the field such as the Peters et al. "invariant prediction" approach.  We'll introduce Pearl's graphical model approach, popular in artificial intelligence and structure learning, in contrast with Rubin's potential outcomes approach, which is more familiar with statisticians.  While the graphical approach emphasizes discovering the entire causal structure of the graph, the potential outcomes approach is more narrowing concerned with estimating a particular causal effect.  The ``invariant prediction'' approach may provide an intermediate approach, which learns part of the graphical structure.

In the first part of the talk, we'll cover the basics of the graphical approach in causality, and illustrate the three main principles.  The causal graph tells you which variables are affected by intervention and also implies certain conditional independence relationships.  Using simple examples, we illustrate the principle of invariant prediction, and derivation of both the back-door criterion and matching approach for estimating an average treatment approach.  Our focus in this first part of the talk is on core principles and theory: we leave practical issues and details (including many criticisms and debates between causal inference researchers, statisticians, and practitioners) for the second half of the talk.
