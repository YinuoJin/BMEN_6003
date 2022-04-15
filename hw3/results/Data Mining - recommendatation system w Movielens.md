#### Data Mining - recommendatation system w/ Movielens



*Collaborative Filtering*

Filter out items that users might like based on similar users' selections.

$\cdot$ smaller set of users with similar preferences $\cdot$ Combine their preferred items $\rightarrow$ ranked list



Model

(1). Memory based

User-based CF: find similar users *having rated* $i$

Item-based CF: find similar items $i$'s that user $u$ *has rated*

Item-based performs more stable than user-based, though performs poorly on entertainment dataset such as Movielens.

- Find similar users/items: based on rating vectors. `scipy.spatial.cosine`<br>

  `similarity = 1 - cosine_dist`

- Factor out user preference: centering ratings 

- Calculate rating $X_{u,i}$:

  - Find user $u$'s top $n$ neighbors $u_{adj}$ with explicit rating on $i$
  - Find the avg. ratings of $u_{adj}$ on $i$: $\hat{X}_{u,i} = \dfrac{\sum_{v \in u_{adj}} X_{v,i}}{\vert u_{adj}\vert}$
  - May give "weighted" score multiplying $X_{v,i} \cdot \text{similar}(u,v)$  

- Accuracy measure: $RMSE(y\_true, y\_pred)$ 



(2). Model based

dimension reduction on user $\times$ item matrix.

Matrix Factorization: $m \times n$ $\rightarrow$ $m \times p$. $p \times n$. $p$: latent feature representation.

- SVD 





*MovieLens* dataset

Matrix form: $X_{u.i}$ denotes user $u$'s rating on item $i$ . Could be sparse, impute missing values given existing ones. 

*MovieLens* summary: 943 users $\times$ 1682 movies

