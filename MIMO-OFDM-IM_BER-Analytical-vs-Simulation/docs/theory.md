# System model and BER analysis

## 1. Transmitter

Each of the $N_t$ transmit antennas sends its own OFDM-IM frame (spatial multiplexing). Each frame has $N_F$ subcarriers, split into $G$ groups of $N = N_F/G$ subcarriers.

In every group:

- $b = \lfloor \log_2 \binom{N}{K} \rfloor$ **index bits** choose one of the first $2^b$ patterns of $K$ active subcarriers (taken from `nchoosek(1:N,K)` in lexicographic order).
- $K\log_2 M$ **symbol bits** are carried by Gray-labelled M-PSK or M-QAM symbols on the active subcarriers.
- In **DM-OFDM-IM**, the $N-K$ idle subcarriers also carry $(N-K)\log_2 M$ bits. They use a second constellation, rotated by $\pi/M$ for PSK and by $\pi/4$ for QAM.

The bits per group are $p = b + K\log_2 M$ for OFDM-IM and $p = b + N\log_2 M$ for DM-OFDM-IM. One group can therefore send $Q = 2^p$ different signals. All of them are stored as the columns of the codebook $\mathbf C \in \mathbb C^{N\times Q}$.

**Interleaving.** Subcarrier $n$ of group $g$ is placed on OFDM subcarrier $(n-1)G + g$, so the subcarriers of one group are $G$ bins apart. This lowers the correlation between their channel gains.

## 2. Channel and receiver

Each link $(r,t)$ is a Rayleigh channel with $L$ taps $h_{r,t}[l] \sim \mathcal{CN}(0, \text{pdp}_l)$, where the power-delay profile satisfies $\sum_l \text{pdp}_l = 1$. Because the cyclic prefix satisfies $CP \ge L-1$, every subcarrier $k$ sees a flat channel:

$$
y_r(k) = \sum_{t=1}^{N_t} H_{r,t}(k)\,x_t(k) + w_r(k),\qquad w_r(k)\sim\mathcal{CN}(0,N_0).
$$

Inside a group, the channel gains are correlated with

$$
\mathbf K_{n_1 n_2} = \sum_{l=0}^{L-1}\text{pdp}_l\, e^{-j2\pi l (n_1-n_2) G / N_F}.
$$

**Joint ML detection** (per group, over all transmit antennas):

$$
(\hat{\mathbf x}_1,\dots,\hat{\mathbf x}_{N_t}) = \arg\min \sum_{r=1}^{N_r}\Big\| \mathbf y_r - \sum_{t=1}^{N_t} \operatorname{diag}(\mathbf h_{r,t})\,\mathbf x_t \Big\|^2 .
$$

## 3. SNR definition

The average energy of one group on one antenna, $E_g$, is computed from the codebook ($E_g = K$ for OFDM-IM with unit-energy symbols).

- `EbN0`: $E_b = E_g (N_F+CP)/(N_F\,p)$, so the energy spent on the CP is counted, and $N_0 = E_b / (E_b/N_0)$.
- `EsN0`: $N_0 = (E_g/N) / (E_s/N_0)$.

## 4. Union bound

The BER is bounded by

$$
P_b \le \frac{1}{N_t\,p\,Q^{N_t}} \sum_{\mathbf X}\sum_{\hat{\mathbf X}} P(\mathbf X\to\hat{\mathbf X})\, e(\mathbf X,\hat{\mathbf X}),
$$

where $e(\cdot,\cdot)$ is the number of bit errors (Hamming distance). It is the sum of the per-antenna Hamming distances.

**Pairwise error probability.** Write $\mathbf E = [\mathbf e_1,\dots,\mathbf e_{N_t}]$ with $\mathbf e_t = \mathbf x_t - \hat{\mathbf x}_t$. Given the channel, the probability is

$$
P(\mathbf X\to\hat{\mathbf X}\mid\mathbf H) = Q\!\left(\sqrt{\tfrac{1}{2N_0}\textstyle\sum_r \|\mathbf\Delta\,\mathbf h_r\|^2}\right),\qquad \mathbf\Delta = [\operatorname{diag}(\mathbf e_1),\dots,\operatorname{diag}(\mathbf e_{N_t})].
$$

Apply $Q(x)\approx \tfrac1{12}e^{-x^2/2}+\tfrac14 e^{-2x^2/3}$ (Chiani et al.). Then average over $\mathbf h_r\sim\mathcal{CN}(\mathbf 0,\mathbf I_{N_t}\otimes\mathbf K)$, which are independent for each receive antenna. With Sylvester's determinant identity this gives

$$
P(\mathbf X\to\hat{\mathbf X}) \approx \frac{1}{12}\det(\mathbf I_N + q_1\mathbf A)^{-N_r} + \frac14\det(\mathbf I_N + q_2\mathbf A)^{-N_r},
$$

$$
\mathbf A = \mathbf\Delta(\mathbf I\otimes\mathbf K)\mathbf\Delta^H = \mathbf K \odot (\mathbf E\mathbf E^H),\qquad q_1=\frac1{4N_0},\; q_2=\frac1{3N_0}.
$$

$\mathbf A$ is always $N\times N$, however many antennas there are.

## 5. Efficient evaluation

1. **Group the single-antenna pairs.** All $Q^2$ pairs $(\mathbf x,\hat{\mathbf x})$ of one antenna are sorted into classes that share the same difference vector $\mathbf e$. A common phase is removed first, because it does not change $\mathbf e\mathbf e^H$. Class $u$ stores how many pairs it holds, $c_u$, and their total Hamming distance, $s_u$.
2. **Combine antennas.** A joint error is a choice of one class per antenna, $(u_1,\dots,u_{N_t})$. It stands for $\prod_t c_{u_t}$ pairs with total bit-error weight
   $w = \sum_t s_{u_t}\prod_{v\ne t} c_{u_v}$. When every antenna is in the zero class, $w = 0$, so no explicit exclusion is needed.
3. **Evaluate all SNRs at once.** $\det(\mathbf I + q\mathbf A) = \sum_{k=0}^N e_k(\mathbf A)\,q^k$. The coefficients $e_k$ are computed once per class combination from $\operatorname{tr}(\mathbf A^k)$ using Newton's identities. Every SNR value then only needs a polynomial evaluation.

If the number of class combinations $U^{N_t}$ is larger than `AN.maxTerms`, the sum is estimated from uniformly sampled combinations (the estimate is unbiased).
