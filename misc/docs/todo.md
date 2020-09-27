# What's still to do

- [ ] Correctness of transform function
    
    - [x] ✅ ```RESOLVED``` It may be incorrect: s. artifacts and index errors. But these seem to decrease with increasing bins.

    - [ ] 🟥 ```TODO``` On Lorenz attractor it works perfectly fine, when it comes to the non-noised. when it's noised the function is misbehaving.

- [x] ✅ ```RESOLVED``` n-dimensional implementation

- [ ] 🟥 ```TODO``` How to solve FPE in python?

![formula](https://render.githubusercontent.com/render/math?math={\frac{\partial}{\partial%20t}p(\vec{x},t+\tau|\vec{x},t)=(-\sum_i\frac{\partial}{\partial%20x_i}D_i^{(1)}(\vec{x},t)+\sum_{ij}\frac{\partial}{\partial%20x_ix_j}D_ij^{(2)}(\vec{x},t))\cdotp(\vec{x},t+\tau|\vec{x},t)})

<!-- $$
    \frac{\partial}{\partial t} p(\vec{x}, t + \tau | \vec{x}, t) = (- \sum_i \frac{\partial}{\partial x_i} D_i^{(1)}(\vec{x}, t) + \sum_{ij} \frac{\partial}{\partial x_i x_j} D_ij^{(2)}(\vec{x}, t)) \cdot p(\vec{x}, t + \tau | \vec{x}, t)
$$ -->

- [ ] 🟥 ```TODO``` How to actually derivate a field? 

![formula](https://render.githubusercontent.com/render/math?math={\frac{\partial}{\partial%20x_i}D_i^{(1)}\text{and}\frac{\partial}{\partial%20x_ix_j}D_ij^{(2)}\text{[1](2)}})

<!-- $$
    \frac{\partial}{\partial x_i} D_i^{(1)} \text{ and } \frac{\partial}{\partial x_i x_j} D_ij^{(2)} \text{ [1](2)}
$$ -->

- [ ] 🟥 ```TODO``` Is it true that 

![formula](https://render.githubusercontent.com/render/math?math={<a%20\cdot%20b>=<a>%20\cdot<b>})

<!-- $$
    < a \cdot b> = < a > \cdot < b >
$$ -->
If it's true, it'll make calculation of $D_{ij}^{(2)}$ easier.