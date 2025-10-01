
#include "arith/R_Q_square.h"
#include "arith/R_Q.h"
#include "arith/Z_Q.h"
#include "util/util.h"
#include "core.h"
#include "endecode.h"


//decomp
template< int LOGQ, int N , int LOG_qdiv>
void decomp(const R_Q_square<  LOGQ, N>& ct, 
	R_Q_square<LOGQ - LOG_qdiv, N>& ct_up, R_Q_square<LOGQ - LOG_qdiv, N>& ct_down) {
	for (int i = 0; i < N; i++) {
		shift_right<LOGQ, LOGQ - LOG_qdiv>(ct[0][i], ct_up[0][i]);
		shift_right<LOGQ, LOGQ - LOG_qdiv>(ct[1][i], ct_up[1][i]);

		Z_Q<LOGQ> temp0, temp1, temp2, temp3;
		
		shift_left<LOGQ - LOG_qdiv, LOGQ>(ct_up[0][i], temp0);
		shift_left<LOGQ - LOG_qdiv, LOGQ>(ct_up[1][i], temp1);

		temp2 = ct[0][i];
		temp2 -= temp0;
		temp3 = ct[1][i];
		temp3 -= temp1;

		resize< LOGQ, LOGQ - LOG_qdiv>(temp2, ct_down[0][i]);
		resize< LOGQ, LOGQ - LOG_qdiv>(temp3, ct_down[1][i]);

	}

}

//decomp_pt
template< int LOGQ, int N, int LOG_qdiv>
void decomp(const R_Q<  LOGQ, N>& pt, R_Q<LOGQ - LOG_qdiv, N>& pt_up, R_Q<LOGQ - LOG_qdiv, N>& pt_down) {
	for (int i = 0; i < N; i++) {
		shift_right<LOGQ, LOGQ - LOG_qdiv>(pt[i], pt_up[i]);

		Z_Q<LOGQ> temp, temp2;

		shift_left<LOGQ - LOG_qdiv, LOGQ>(pt_up[i], temp);

		temp2 = pt[i];
		temp2 -= temp;

		resize< LOGQ, LOGQ - LOG_qdiv>(temp2, pt_down[i]);

	}

}

//recomb
template< int LOGQ, int N, int LOG_qdiv>
void recomb(const R_Q_square<LOGQ , N>& ct_up, const R_Q_square<LOGQ, N>& ct_down
	, R_Q_square<  LOGQ, N>& ct_comb) {
	for (int i = 0; i < N; i++) {

		Z_Q<LOGQ > temp0, temp1;
		temp0 = ct_up[0][i];
		temp0 *= 1 << LOG_qdiv;
		temp0 += ct_down[0][i];
		ct_comb[0][i] = temp0;

		temp1 = ct_up[1][i];
		temp1 *= 1 << LOG_qdiv;
		temp1 += ct_down[1][i];
		ct_comb[1][i] = temp1;
	}

}

//tuple_add
template< int LOGQ, int N>
void tuple_add(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, const R_Q_square<LOGQ, N>& ct_up2, const R_Q_square<LOGQ, N>& ct_down2
    , R_Q_square<LOGQ, N>& ct_up3, R_Q_square<LOGQ, N>& ct_down3) {
	for (int i = 0; i < N; i++) {
		ct_up3[0][i] = ct_up1[0][i];
		ct_up3[0][i] += ct_up2[0][i];
		ct_up3[1][i] = ct_up1[1][i];
		ct_up3[1][i] += ct_up2[1][i];

		ct_down3[0][i] = ct_down1[0][i];
		ct_down3[0][i] += ct_down2[0][i];
		ct_down3[1][i] = ct_down1[1][i];
		ct_down3[1][i] += ct_down2[1][i];
	}

}

//tuple_tensor_c
template< int LOGQ, int N>
void tuple_tensor(const R_Q<LOGQ , N>& pt_up, const R_Q<LOGQ , N>& pt_down,
	const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, R_Q_square<LOGQ, N>& ct_up2, R_Q_square<LOGQ, N>& ct_down2) {
	
	conv<LOGQ, N>(pt_up, ct_up1[0], ct_up2[0]);
	conv<LOGQ, N>(pt_up, ct_up1[1], ct_up2[1]);

	R_Q<LOGQ, N> temp;
	conv<LOGQ, N>(pt_down, ct_up1[0], temp);
	ct_down2[0] = temp;
	conv<LOGQ, N>(pt_up, ct_down1[0], temp);
	ct_down2[0] += temp;

	conv<LOGQ, N>(pt_down, ct_up1[1], temp);
	ct_down2[1] = temp;
	conv<LOGQ, N>(pt_up, ct_down1[1], temp);
	ct_down2[1] += temp;

}

//tuple_rs
template< int LOGQ, int N, int LOG_qdiv, int LOG_ql>
void tuple_rs(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, R_Q_square<LOGQ - LOG_ql, N>& ct_up2, R_Q_square<LOGQ - LOG_ql, N>& ct_down2) {

	RS<LOGQ, LOGQ - LOG_ql, N>(ct_up1, ct_up2);
	
	R_Q_square<LOGQ, N> ct_comb;
	recomb<LOGQ, N, LOG_qdiv>(ct_up1, ct_down1, ct_comb);

	RS<LOGQ, LOGQ - LOG_ql, N>(ct_comb, ct_down2);

	R_Q_square<LOGQ - LOG_ql, N> temp;
	RS<LOGQ, LOGQ - LOG_ql, N>(ct_up1, temp);

	temp *= 1 << LOG_qdiv;

	ct_down2 -= temp;

}

//tuple_cmult
template< int LOGQ, int N, int LOG_qdiv, int LOG_ql>
void tuple_cmult(const R_Q<LOGQ, N>& pt_up, const R_Q<LOGQ, N>& pt_down,
	const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, R_Q_square<LOGQ - LOG_ql, N>& ct_up2, R_Q_square<LOGQ - LOG_ql, N>& ct_down2) {

	R_Q_square<LOGQ, N> temp_up, temp_down;
	tuple_tensor<LOGQ, N>(pt_up, pt_down, ct_up1, ct_down1, temp_up, temp_down);

	tuple_rs< LOGQ, N, LOG_qdiv, LOG_ql>(temp_up, temp_down, ct_up2, ct_down2);

}


//tuple_permute
template< int LOGQ, int N, int LOG_qdiv>
void tuple_permute(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1, const int j
	, R_Q_square<LOGQ, N>& ct_up2, R_Q_square<LOGQ, N>& ct_down2) {

	rot<LOGQ, N>(ct_up1[0], ct_up2[0], j);
	rot<LOGQ, N>(ct_up1[1], ct_up2[1], j);
	rot<LOGQ, N>(ct_down1[0], ct_down2[0], j);
	rot<LOGQ, N>(ct_down1[1], ct_down2[1], j);

}

//tuple_conj
template< int LOGQ, int N, int LOG_qdiv>
void tuple_conj_pre(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, R_Q_square<LOGQ, N>& ct_up2, R_Q_square<LOGQ, N>& ct_down2) {

	conj<LOGQ, N>(ct_up1[0], ct_up2[0]);
	conj<LOGQ, N>(ct_up1[1], ct_up2[1]);
	conj<LOGQ, N>(ct_down1[0], ct_down2[0]);
	conj<LOGQ, N>(ct_down1[1], ct_down2[1]);

}


//tuple_ksw
template< int LOGQ, int N, int LOG_qdiv>
void tuple_ksw(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, const R_Q_square<2 * (LOGQ+ LOG_qdiv), N>& rkey1, const R_Q_square<2 * LOGQ , N>& rkey2, 
	R_Q_square<LOGQ, N>& ct_up2, R_Q_square<LOGQ, N>& ct_down2) {

	R_Q_square<LOGQ + LOG_qdiv, N> qdiv_ct_up;
	Inv_RS<LOGQ, LOGQ + LOG_qdiv, N>(ct_up1, qdiv_ct_up);

	HEAAN<LOGQ + LOG_qdiv, N>::ks(rkey1, qdiv_ct_up);
	decomp<LOGQ + LOG_qdiv, N, LOG_qdiv>(qdiv_ct_up, ct_up2, ct_down2);

	R_Q_square<LOGQ , N> temp;
	temp = ct_down1;
	HEAAN<LOGQ, N>::ks(rkey2, temp);
	ct_down2 += temp;

}



//tuple_rot
template< int LOGQ, int N, int LOG_qdiv>
void tuple_rot(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1, const int j
	, const R_Q_square<2 * (LOGQ + LOG_qdiv), N>& rkey1, const R_Q_square<2 * LOGQ, N>& rkey2
	, R_Q_square<LOGQ, N>& ct_up2, R_Q_square<LOGQ, N>& ct_down2) {

	R_Q_square<LOGQ, N> temp_up, temp_down;
	tuple_permute<LOGQ, N, LOG_qdiv>(ct_up1, ct_down1, j, temp_up, temp_down);
	tuple_ksw<LOGQ, N, LOG_qdiv>(temp_up, temp_down, rkey1, rkey2, ct_up2, ct_down2);

}

//tuple_conj
template< int LOGQ, int N, int LOG_qdiv>
void tuple_conj_full(const R_Q_square<LOGQ, N>& ct_up1, const R_Q_square<LOGQ, N>& ct_down1
	, const R_Q_square<2 * (LOGQ + LOG_qdiv), N>& rkey1, const R_Q_square<2 * LOGQ, N>& rkey2
	, R_Q_square<LOGQ, N>& ct_up2, R_Q_square<LOGQ, N>& ct_down2) {

	R_Q_square<LOGQ, N> temp_up, temp_down;
	tuple_conj_pre<LOGQ, N, LOG_qdiv>(ct_up1, ct_down1, temp_up, temp_down);
	tuple_ksw<LOGQ, N, LOG_qdiv>(temp_up, temp_down, rkey1, rkey2, ct_up2, ct_down2);

}

//tuple_lintrans
template< int LOGQ, int LOGN, int LOGDELTA, int S, int LOG_qdiv >
void tuple_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S>& Ar,
	const SparseDiagonal<1 << (LOGN - 1), S>& Ai,
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	const int N = 1 << LOGN;
	Act_up.setzero();
	Act_down.setzero();
	for (int s = 0; s < S; s++) {
		assert(Ar.off[s] == Ai.off[s]);
		if (Ar.zero[s] && Ai.zero[s])
			continue;
		R_Q<LOGQ + LOG_qdiv, N> pt;
		encode<LOGQ + LOG_qdiv, LOGN>(Ar.vec[s],
			Ai.vec[s], 1ULL << LOGDELTA, pt);

		R_Q<LOGQ, N> pt_up, pt_down;
		decomp<LOGQ + LOG_qdiv, N, LOG_qdiv>(pt, pt_up, pt_down);

		R_Q_square<LOGQ, N> ct_rot_up(ct_up), ct_rot_down(ct_down);
		if (Ar.off[s] != 0) {
			R_Q_square<2 * (LOGQ+LOG_qdiv), 1 << LOGN> rkey1;
			R_Q_square<2 * LOGQ, 1 << LOGN> rkey2;
			int skey_rot[N];
			rot<N>(skey, skey_rot, Ar.off[s]);
			HEAAN<LOGQ + LOG_qdiv, N>::swkgen(skey_rot, skey, rkey1);
			HEAAN<LOGQ, N>::swkgen(skey_rot, skey, rkey2);

			tuple_rot<LOGQ, N, LOG_qdiv>(ct_up, ct_down, Ar.off[s], rkey1, rkey2, ct_rot_up, ct_rot_down);
		}

		R_Q_square<LOGQ, N> tensor_up, tensor_down;
		tuple_tensor<LOGQ, N>(pt_up, pt_down, ct_rot_up, ct_rot_down, tensor_up, tensor_down);

		R_Q_square<LOGQ, N> temp_up, temp_down;
		temp_up = Act_up;
		temp_down = Act_down;
		
		tuple_add(temp_up, temp_down, tensor_up, tensor_down, Act_up, Act_down);
	}
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	R_Q_square<LOGQ, 1 << LOGN> ct_temp_up, ct_temp_down;
	Act_up = ct_up;
	Act_down = ct_down;
	for (int d = 0; d < D; ++d) {
		ct_temp_up = Act_up;
		ct_temp_down = Act_down;
		tuple_linear_transform<LOGQ, LOGN, LOGDELTA, S, LOG_qdiv>(Ar[d], Ai[d], ct_temp_up, ct_temp_down, skey, 
			Act_up, Act_down);
	}
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped2_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	assert(D % 2 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

	for (int i = 0; i < D / 2; ++i) {
		MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
	}
	tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S, D / 2, LOG_qdiv>(Br, Bi, ct_up, ct_down, skey, Act_up, Act_down);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped3_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	assert(D % 3 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 3];

	for (int i = 0; i < D / 3; ++i) {
		MatMul(Ar[3 * i + 1], Ai[3 * i + 1], Ar[3 * i], Ai[3 * i], Br[i], Bi[i]);
		MatMul(Ar[3 * i + 2], Ai[3 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
	}
	tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S, D / 3, LOG_qdiv>(Cr, Ci, ct_up, ct_down, skey, Act_up, Act_down);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped4_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	if(D % 4 == 0){
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 4];

	for (int i = 0; i < D / 4; ++i) {
		MatMul(Ar[4 * i + 1], Ai[4 * i + 1], Ar[4 * i], Ai[4 * i], Br[i], Bi[i]);
		MatMul(Ar[4 * i + 2], Ai[4 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[4 * i + 3], Ai[4 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
	}
	tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4, LOG_qdiv>(Dr, Di, ct_up, ct_down, skey, Act_up, Act_down);
}

	if(D %4== 3){
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 4];

	SparseDiagonal<1 << (LOGN - 1), S* S> LBr;
	SparseDiagonal<1 << (LOGN - 1), S* S> LBi;
	SparseDiagonal<1 << (LOGN - 1), S* S* S> LCr;
	SparseDiagonal<1 << (LOGN - 1), S* S* S> LCi;

	for (int i = 0; i < D / 4; ++i) {
		MatMul(Ar[4 * i + 1], Ai[4 * i + 1], Ar[4 * i], Ai[4 * i], Br[i], Bi[i]);
		MatMul(Ar[4 * i + 2], Ai[4 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[4 * i + 3], Ai[4 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
	}

	MatMul(Ar[D-2], Ai[D-2], Ar[D-3], Ai[D-3], LBr, LBi);
	MatMul(Ar[D-1], Ai[D-1], LBr, LBi, LCr, LCi);

	R_Q_square<LOGQ, 1 << LOGN> temp_up, temp_down;
	tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4, LOG_qdiv>(Dr, Di, ct_up, ct_down, skey, temp_up, temp_down);
	tuple_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S, LOG_qdiv>(LCr, LCi, temp_up, temp_down, skey, Act_up, Act_down);
}

}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped5_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	assert(D % 5 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Er[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Ei[D / 5];

	for (int i = 0; i < D / 5; ++i) {
		MatMul(Ar[5 * i + 1], Ai[5 * i + 1], Ar[5 * i], Ai[5 * i], Br[i], Bi[i]);
		MatMul(Ar[5 * i + 2], Ai[5 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[5 * i + 3], Ai[5 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
		MatMul(Ar[5 * i + 4], Ai[5 * i + 4], Dr[i], Di[i], Er[i], Ei[i]);
	}

	tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D / 5, LOG_qdiv>(Er, Ei, ct_up, ct_down, skey, Act_up, Act_down);

}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped7_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	assert(D % 7 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Er[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Ei[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S> Fr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S> Fi[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S* S> Gr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S* S> Gi[D / 7];
	for (int i = 0; i < D / 7; ++i) {
		MatMul(Ar[7 * i + 1], Ai[7 * i + 1], Ar[7 * i], Ai[7 * i], Br[i], Bi[i]);
		MatMul(Ar[7 * i + 2], Ai[7 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[7 * i + 3], Ai[7 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
		MatMul(Ar[7 * i + 4], Ai[7 * i + 4], Dr[i], Di[i], Er[i], Ei[i]);
		MatMul(Ar[7 * i + 5], Ai[7 * i + 5], Er[i], Ei[i], Fr[i], Fi[i]);
		MatMul(Ar[7 * i + 6], Ai[7 * i + 6], Fr[i], Fi[i], Gr[i], Gi[i]);
	}
	tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D / 7, LOG_qdiv>(Gr, Gi, ct_up, ct_down, skey, Act_up, Act_down);
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G = 1, int LOG_qdiv>
void tuple_grouped_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act_up, R_Q_square<  LOGQ, 1 << LOGN>& Act_down) {
	if (G == 1)
		tuple_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 2)
		tuple_grouped2_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 3)
		tuple_grouped3_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 4)
		tuple_grouped4_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 5)
		tuple_grouped5_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 7)
		tuple_grouped7_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else
		assert(false);
}


template< int LOGQ, int LOGN, int LOGDELTA, int S, int LOG_qdiv, int LOG_ql>
void tuple_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S>& Ar,
	const SparseDiagonal<1 << (LOGN - 1), S>& Ai,
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	const int N = 1 << LOGN;
	Act_up.setzero(); Act_down.setzero();
	for (int s = 0; s < S; s++) {
		assert(Ar.off[s] == Ai.off[s]);
		if (Ar.zero[s] && Ai.zero[s])
			continue;
		R_Q<LOGQ + LOG_qdiv, N> pt;
		encode<LOGQ + LOG_qdiv, LOGN>(Ar.vec[s],
			Ai.vec[s], 1ULL << LOGDELTA, pt);

		R_Q<LOGQ + LOG_qdiv, N> pt_rot;
		R_Q<LOGQ, N> pt_up, pt_down;

		R_Q_square<LOGQ, N> ct_rot_up, ct_rot_down;
		R_Q_square<LOGQ - LOG_ql, N> ct_block_up, ct_block_down;

		if (Ar.off[s] != 0) {

			rot<LOGQ + LOG_qdiv, N>(pt, pt_rot, N / 2 - Ar.off[s]);
			decomp<LOGQ + LOG_qdiv, N, LOG_qdiv>(pt_rot, pt_up, pt_down);
			tuple_tensor<LOGQ, N>(pt_up, pt_down, ct_up, ct_down, ct_rot_up, ct_rot_down);

			R_Q_square<LOGQ - LOG_ql, N> ct_rs_up, ct_rs_down;

			tuple_rs<LOGQ, N, LOG_qdiv, LOG_ql>(ct_rot_up, ct_rot_down, ct_rs_up, ct_rs_down);

			R_Q_square<2 * (LOGQ + LOG_qdiv - LOG_ql), 1 << LOGN> rkey1;
			R_Q_square<2 * (LOGQ - LOG_ql), 1 << LOGN> rkey2;
			int skey_rot[N];
			rot<N>(skey, skey_rot, Ar.off[s]);
			HEAAN<LOGQ + LOG_qdiv - LOG_ql, N>::swkgen(skey_rot, skey, rkey1);
			HEAAN<LOGQ - LOG_ql, N>::swkgen(skey_rot, skey, rkey2);

			tuple_rot<LOGQ - LOG_ql, N, LOG_qdiv>(ct_rs_up, ct_rs_down, Ar.off[s], rkey1, rkey2, ct_block_up, ct_block_down);
		}

		if (Ar.off[s] == 0) {
			decomp<LOGQ + LOG_qdiv, N, LOG_qdiv>(pt, pt_up, pt_down);
			tuple_tensor<LOGQ, N>(pt_up, pt_down, ct_up, ct_down, ct_rot_up, ct_rot_down);
			tuple_rs<LOGQ, N, LOG_qdiv, LOG_ql>(ct_rot_up, ct_rot_down, ct_block_up, ct_block_down);
		}

		tuple_add(Act_up, Act_down, ct_block_up, ct_block_down, Act_up, Act_down);
	}
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv, int LOG_ql>
void tuple_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {

	Act_up.setzero(); Act_down.setzero();

	tuple_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, LOG_qdiv, LOG_ql>(Ar[0], Ai[0], ct_up, ct_down, skey, Act_up, Act_down);

	R_Q_square<LOGQ - LOG_ql, 1 << LOGN> ct_temp_up, ct_temp_down;

	for (int d = 1; d < D; ++d) {
		ct_temp_up = Act_up;
		ct_temp_down = Act_down;
		tuple_linear_transform<LOGQ - LOG_ql, LOGN, LOGDELTA, S, LOG_qdiv>(Ar[d], Ai[d], ct_temp_up, ct_temp_down, skey, Act_up, Act_down);
	}
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv, int LOG_ql>
void tuple_grouped2_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	assert(D % 2 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

	for (int i = 0; i < D / 2; ++i) {
		MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
	}
	tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S, D / 2, LOG_qdiv, LOG_ql>(Br, Bi, ct_up, ct_down, skey, Act_up, Act_down);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv, int LOG_ql>
void tuple_grouped3_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	assert(D % 3 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 3];

	for (int i = 0; i < D / 3; ++i) {
		MatMul(Ar[3 * i + 1], Ai[3 * i + 1], Ar[3 * i], Ai[3 * i], Br[i], Bi[i]);
		MatMul(Ar[3 * i + 2], Ai[3 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
	}
	tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S, D / 3, LOG_qdiv, LOG_ql>(Cr, Ci, ct_up, ct_down, skey, Act_up, Act_down);
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv, int LOG_ql>
void tuple_grouped4_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	if(D % 4 == 0){
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 4];

	for (int i = 0; i < D / 4; ++i) {
		MatMul(Ar[4 * i + 1], Ai[4 * i + 1], Ar[4 * i], Ai[4 * i], Br[i], Bi[i]);
		MatMul(Ar[4 * i + 2], Ai[4 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[4 * i + 3], Ai[4 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
	}
	tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4, LOG_qdiv, LOG_ql>(Dr, Di, ct_up, ct_down, skey, Act_up, Act_down);
}

	if(D % 4 == 3){
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 4];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 4];

	SparseDiagonal<1 << (LOGN - 1), S* S> LBr;
	SparseDiagonal<1 << (LOGN - 1), S* S> LBi;
	SparseDiagonal<1 << (LOGN - 1), S* S* S> LCr;
	SparseDiagonal<1 << (LOGN - 1), S* S* S> LCi;

	for (int i = 0; i < D / 4; ++i) {
		MatMul(Ar[4 * i + 1], Ai[4 * i + 1], Ar[4 * i], Ai[4 * i], Br[i], Bi[i]);
		MatMul(Ar[4 * i + 2], Ai[4 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[4 * i + 3], Ai[4 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
	}

	MatMul(Ar[D-2], Ai[D-2], Ar[D-3], Ai[D-3], LBr, LBi);
	MatMul(Ar[D-1], Ai[D-1], LBr, LBi, LCr, LCi);

	R_Q_square<LOGQ - LOG_ql, 1 << LOGN> temp_up, temp_down;
	tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4, LOG_qdiv, LOG_ql>(Dr, Di, ct_up, ct_down, skey, temp_up, temp_down);
	tuple_linear_transform<LOGQ - LOG_ql, LOGN, LOGDELTA, S* S* S, LOG_qdiv>(LCr, LCi, temp_up, temp_down, skey, Act_up, Act_down);
}

}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv, int LOG_ql>
void tuple_grouped5_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	assert(D % 5 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Er[D / 5];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Ei[D / 5];

	for (int i = 0; i < D / 5; ++i) {
		MatMul(Ar[5 * i + 1], Ai[5 * i + 1], Ar[5 * i], Ai[5 * i], Br[i], Bi[i]);
		MatMul(Ar[5 * i + 2], Ai[5 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[5 * i + 3], Ai[5 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
		MatMul(Ar[5 * i + 4], Ai[5 * i + 4], Dr[i], Di[i], Er[i], Ei[i]);
	}
	tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D / 5, LOG_qdiv, LOG_ql>(Er, Ei, ct_up, ct_down, skey, Act_up, Act_down);
}


template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv, int LOG_ql>
void tuple_grouped7_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	assert(D % 7 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Er[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S> Ei[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S> Fr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S> Fi[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S* S> Gr[D / 7];
	SparseDiagonal<1 << (LOGN - 1), S* S* S* S* S* S* S> Gi[D / 7];
	for (int i = 0; i < D / 7; ++i) {
		MatMul(Ar[7 * i + 1], Ai[7 * i + 1], Ar[7 * i], Ai[7 * i], Br[i], Bi[i]);
		MatMul(Ar[7 * i + 2], Ai[7 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
		MatMul(Ar[7 * i + 3], Ai[7 * i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
		MatMul(Ar[7 * i + 4], Ai[7 * i + 4], Dr[i], Di[i], Er[i], Ei[i]);
		MatMul(Ar[7 * i + 5], Ai[7 * i + 5], Er[i], Ei[i], Fr[i], Fi[i]);
		MatMul(Ar[7 * i + 6], Ai[7 * i + 6], Fr[i], Fi[i], Gr[i], Gi[i]);
	}
	tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D / 7, LOG_qdiv, LOG_ql>(Gr, Gi, ct_up, ct_down, skey, Act_up, Act_down);
}


template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G, int LOG_qdiv, int LOG_ql >
void tuple_grouped_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct_up, const R_Q_square<  LOGQ, 1 << LOGN>& ct_down,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_up, R_Q_square<  LOGQ - LOG_ql, 1 << LOGN>& Act_down) {
	if (G == 1)
		tuple_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv, LOG_ql>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 2)
		tuple_grouped2_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv, LOG_ql>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 3)
		tuple_grouped3_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv, LOG_ql>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 4)
		tuple_grouped4_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv, LOG_ql>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 5)
		tuple_grouped5_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv, LOG_ql>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else if (G == 7)
		tuple_grouped7_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv, LOG_ql>(Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down);
	else
		assert(false);
}





















// ============================================================================
// Tuple Linear Transform with Baby-Step Giant-Step (BSGS) + Adaptive wrappers
// ============================================================================

template<int LOGQ, int LOGN, int LOGDELTA, int S, int LOG_qdiv>
void tuple_linear_transform_sparse_babystep_giantstep(
    const SparseDiagonal<(1<<(LOGN-1)), S>& Ar,
    const SparseDiagonal<(1<<(LOGN-1)), S>& Ai,
    const R_Q_square<LOGQ, (1<<LOGN)>& ct_up,
    const R_Q_square<LOGQ, (1<<LOGN)>& ct_down,
    const int skey[(1<<LOGN)],
    R_Q_square<LOGQ, (1<<LOGN)>& Act_up,
    R_Q_square<LOGQ, (1<<LOGN)>& Act_down,
    int baby_step_size 
){
    constexpr int N    = 1 << LOGN;
    constexpr int HALF = N >> 1;

    auto normH = [](int x)->int{
        int y = x % HALF;
        return (y < 0) ? (y + HALF) : y;
    };
    auto igcd = [](int a, int b){
        a = std::abs(a); b = std::abs(b);
        while (b){ int t = a % b; a = b; b = t; }
        return a;
    };
    auto get_off = [&](int s)->int {
        if (!Ar.zero[s] && !Ai.zero[s]){
            assert(Ar.off[s] == Ai.off[s] && "Ar.off[s] != Ai.off[s]");
            return Ar.off[s];
        } else if (!Ar.zero[s]) {
            return Ar.off[s];
        } else if (!Ai.zero[s]) {
            return Ai.off[s];
        } else {
            return 0;
        }
    };

    std::vector<int> offs; offs.reserve(S);
    Act_up.setzero();
    Act_down.setzero();
    for (int s=0; s<S; ++s){
        if (!Ar.zero[s] || !Ai.zero[s]){
            offs.push_back(normH(get_off(s)));
        }
    }
    if (offs.empty()) return;

    int k = 0;
    for (int v: offs) if (v) k = (k==0 ? v : igcd(k, v));
    if (k == 0) k = 1;
    std::cerr << "[Tuple BSGS] stride k = " << k << '\n';

    for (int v: offs) { assert(normH(v) % k == 0 && "offset not multiple of stride k"); }

    std::vector<int> opOf(S, 0);
    std::set<int> uniqOp;
    for (int s=0; s<S; ++s){
        if (Ar.zero[s] && Ai.zero[s]) continue;
        int op = normH(get_off(s)) / k;
        opOf[s] = op;
        uniqOp.insert(op);
    }

    int M  = (int)uniqOp.size();
    int N1 = (baby_step_size > 0 ? baby_step_size
                                 : std::max(1, (int)std::sqrt((double)M)));

    std::map<int, std::vector<int>> groupJ; 
    std::set<int> babyI;
    for (int s=0; s<S; ++s){
        if (Ar.zero[s] && Ai.zero[s]) continue;
        int op = opOf[s];
        int j  = op / N1;
        int i  = op % N1;
        groupJ[j].push_back(s);
        babyI.insert(i);
    }

    std::map<int, R_Q_square<LOGQ, N>> ct_baby_up;
    std::map<int, R_Q_square<LOGQ, N>> ct_baby_down;
    for (int i : babyI){
        int amt = normH(i * k);
        if (amt == 0){
            ct_baby_up[i]   = ct_up;
            ct_baby_down[i] = ct_down;
            continue;
        }
        R_Q_square<2*(LOGQ+LOG_qdiv), N> rkey1;
        R_Q_square<2*LOGQ, N>         rkey2;
        int skey_rot[N];
        rot<N>(skey, skey_rot, amt);
        HEAAN<LOGQ+LOG_qdiv, N>::swkgen(skey_rot, skey, rkey1);
        HEAAN<LOGQ, N>::swkgen(skey_rot, skey, rkey2);

        R_Q_square<LOGQ, N> up_rot, down_rot;
        tuple_rot<LOGQ, N, LOG_qdiv>(ct_up, ct_down, amt, rkey1, rkey2, up_rot, down_rot);
        ct_baby_up[i]   = up_rot;
        ct_baby_down[i] = down_rot;
    }
    if (!ct_baby_up.count(0)){ ct_baby_up[0] = ct_up; }
    if (!ct_baby_down.count(0)){ ct_baby_down[0] = ct_down; }

    for (auto &kv : groupJ){
        const int j = kv.first;
        const auto& idxs = kv.second;

        R_Q_square<LOGQ, N> inner_up;   inner_up.setzero();
        R_Q_square<LOGQ, N> inner_down; inner_down.setzero();

        const int r = normH(N1 * j * k);

        for (int sIdx : idxs){
            int i = opOf[sIdx] % N1;

            std::vector<double> re(HALF,0.0), im(HALF,0.0);
            if (!Ar.zero[sIdx]){
                const double* src = Ar.vec[sIdx];
                for (int t=0; t<HALF; ++t) re[t] = src[(t + HALF - r) % HALF];
            }
            if (!Ai.zero[sIdx]){
                const double* src = Ai.vec[sIdx];
                for (int t=0; t<HALF; ++t) im[t] = src[(t + HALF - r) % HALF];
            }

            R_Q<LOGQ + LOG_qdiv, N> pt_all;
            encode<LOGQ + LOG_qdiv, LOGN>(re.data(), im.data(), 1ULL<<LOGDELTA, pt_all);

            R_Q<LOGQ, N> pt_up, pt_down;
            decomp<LOGQ + LOG_qdiv, N, LOG_qdiv>(pt_all, pt_up, pt_down);

            R_Q_square<LOGQ, N> tens_up, tens_down;
            const auto& baby_up   = ct_baby_up[i];
            const auto& baby_down = ct_baby_down[i];
            tuple_tensor<LOGQ, N>(pt_up, pt_down, baby_up, baby_down, tens_up, tens_down);

            R_Q_square<LOGQ, N> tmp_up = inner_up, tmp_down = inner_down;
            tuple_add(tmp_up, tmp_down, tens_up, tens_down, inner_up, inner_down);
        }

        if (j != 0){
            int amt = r; 
            R_Q_square<2*(LOGQ+LOG_qdiv), N> rkey1;
            R_Q_square<2*LOGQ, N>         rkey2;
            int skey_rot[N];
            rot<N>(skey, skey_rot, amt);
            HEAAN<LOGQ+LOG_qdiv, N>::swkgen(skey_rot, skey, rkey1);
            HEAAN<LOGQ, N>::swkgen(skey_rot, skey, rkey2);

            R_Q_square<LOGQ, N> inner_rot_up, inner_rot_down;
            tuple_rot<LOGQ, N, LOG_qdiv>(inner_up, inner_down, amt, rkey1, rkey2,
                                         inner_rot_up, inner_rot_down);

            R_Q_square<LOGQ, N> tmp_up = Act_up, tmp_down = Act_down;
            tuple_add(tmp_up, tmp_down, inner_rot_up, inner_rot_down, Act_up, Act_down);
        } else {
            R_Q_square<LOGQ, N> tmp_up = Act_up, tmp_down = Act_down;
            tuple_add(tmp_up, tmp_down, inner_up, inner_down, Act_up, Act_down);
        }
    }
}


// ----------------------------------------------------------------------------
// Adaptive wrapper
// ----------------------------------------------------------------------------
template<int LOGQ, int LOGN, int LOGDELTA, int S, int LOG_qdiv>
void tuple_linear_transform_adaptive_sparse_babystep_giantstep(
    const SparseDiagonal<1<<(LOGN-1), S>& Ar,
    const SparseDiagonal<1<<(LOGN-1), S>& Ai,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    int cnt = 0;
    for (int s=0; s<S; ++s) if (!Ar.zero[s] || !Ai.zero[s]) ++cnt;
    int suggested_baby = std::max(1, (int)std::sqrt((double)cnt));

    tuple_linear_transform_sparse_babystep_giantstep<LOGQ, LOGN, LOGDELTA, S, LOG_qdiv>(
        Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down, suggested_baby
    );
}


// ----------------------------------------------------------------------------
// Serial application of D transforms (adaptive BSGS inside)
// ----------------------------------------------------------------------------
template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    R_Q_square<LOGQ, 1<<LOGN> ct_tmp_up, ct_tmp_down;
    Act_up   = ct_up;
    Act_down = ct_down;
    for (int d=0; d<D; ++d){
        ct_tmp_up   = Act_up;
        ct_tmp_down = Act_down;
        tuple_linear_transform_adaptive_sparse_babystep_giantstep<LOGQ, LOGN, LOGDELTA, S, LOG_qdiv>(
            Ar[d], Ai[d], ct_tmp_up, ct_tmp_down, skey, Act_up, Act_down
        );
    }
}


// ----------------------------------------------------------------------------
// Grouped G = 2,3,4(0/3),5,7  (adaptive BSGS inside)
// ----------------------------------------------------------------------------
template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped2_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    assert(D % 2 == 0);
    SparseDiagonal<1<<(LOGN-1), S* S> Br[D/2];
    SparseDiagonal<1<<(LOGN-1), S* S> Bi[D/2];
    for (int i=0; i<D/2; ++i){
        MatMul(Ar[2*i+1], Ai[2*i+1], Ar[2*i], Ai[2*i], Br[i], Bi[i]);
    }
    tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S, D/2, LOG_qdiv>(
        Br, Bi, ct_up, ct_down, skey, Act_up, Act_down
    );
}


template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped3_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    assert(D % 3 == 0);
    SparseDiagonal<1<<(LOGN-1), S* S> Br[D/3];
    SparseDiagonal<1<<(LOGN-1), S* S> Bi[D/3];
    SparseDiagonal<1<<(LOGN-1), S* S* S> Cr[D/3];
    SparseDiagonal<1<<(LOGN-1), S* S* S> Ci[D/3];
    for (int i=0; i<D/3; ++i){
        MatMul(Ar[3*i+1], Ai[3*i+1], Ar[3*i], Ai[3*i], Br[i], Bi[i]);
        MatMul(Ar[3*i+2], Ai[3*i+2], Br[i], Bi[i], Cr[i], Ci[i]);
    }
    tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S, D/3, LOG_qdiv>(
        Cr, Ci, ct_up, ct_down, skey, Act_up, Act_down
    );
}


template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped4_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    if (D == 0){ Act_up.setzero(); Act_down.setzero(); return; }

    if (D % 4 == 0){
        SparseDiagonal<1<<(LOGN-1), S* S>       Br[D/4];
        SparseDiagonal<1<<(LOGN-1), S* S>       Bi[D/4];
        SparseDiagonal<1<<(LOGN-1), S* S* S>    Cr[D/4];
        SparseDiagonal<1<<(LOGN-1), S* S* S>    Ci[D/4];
        SparseDiagonal<1<<(LOGN-1), S* S* S* S> Dr[D/4];
        SparseDiagonal<1<<(LOGN-1), S* S* S* S> Di[D/4];

        for (int i=0; i<D/4; ++i){
            MatMul(Ar[4*i+1], Ai[4*i+1], Ar[4*i],     Ai[4*i],     Br[i], Bi[i]);
            MatMul(Ar[4*i+2], Ai[4*i+2], Br[i],       Bi[i],       Cr[i], Ci[i]);
            MatMul(Ar[4*i+3], Ai[4*i+3], Cr[i],       Ci[i],       Dr[i], Di[i]);
        }
        tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S, D/4, LOG_qdiv>(
            Dr, Di, ct_up, ct_down, skey, Act_up, Act_down
        );
        return;
    }

    if (D % 4 == 3){
        SparseDiagonal<1<<(LOGN-1), S* S>    LBr;
        SparseDiagonal<1<<(LOGN-1), S* S>    LBi;
        SparseDiagonal<1<<(LOGN-1), S* S* S> LCr;
        SparseDiagonal<1<<(LOGN-1), S* S* S> LCi;

        MatMul(Ar[1], Ai[1], Ar[0], Ai[0], LBr, LBi);
        MatMul(Ar[2], Ai[2], LBr,  LBi,    LCr, LCi);

        R_Q_square<LOGQ, 1<<LOGN> ct_after3_up, ct_after3_down;
        tuple_linear_transform_adaptive_sparse_babystep_giantstep<LOGQ, LOGN, LOGDELTA, S* S* S, LOG_qdiv>(
            LCr, LCi, ct_up, ct_down, skey, ct_after3_up, ct_after3_down
        );

        if (D == 3){ Act_up = ct_after3_up; Act_down = ct_after3_down; return; }

        constexpr int R = (D - 3) / 4;
        SparseDiagonal<1<<(LOGN-1), S* S>       Br[R];
        SparseDiagonal<1<<(LOGN-1), S* S>       Bi[R];
        SparseDiagonal<1<<(LOGN-1), S* S* S>    Cr[R];
        SparseDiagonal<1<<(LOGN-1), S* S* S>    Ci[R];
        SparseDiagonal<1<<(LOGN-1), S* S* S* S> Dr[R];
        SparseDiagonal<1<<(LOGN-1), S* S* S* S> Di[R];

        for (int i=0; i<R; ++i){
            const int base = 3 + 4*i;
            MatMul(Ar[base+1], Ai[base+1], Ar[base+0], Ai[base+0], Br[i], Bi[i]);
            MatMul(Ar[base+2], Ai[base+2], Br[i],      Bi[i],      Cr[i], Ci[i]);
            MatMul(Ar[base+3], Ai[base+3], Cr[i],      Ci[i],      Dr[i], Di[i]);
        }
        tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S, R, LOG_qdiv>(
            Dr, Di, ct_after3_up, ct_after3_down, skey, Act_up, Act_down
        );
        return;
    }

    assert(false && "tuple_grouped4_serial_linear_transform_adaptive: D % 4 must be 0 or 3");
}


template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped5_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    assert(D % 5 == 0);
    SparseDiagonal<1<<(LOGN-1), S* S>             Br[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S>             Bi[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S* S>          Cr[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S* S>          Ci[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S>       Dr[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S>       Di[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S>    Er[D/5];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S>    Ei[D/5];

    for (int i=0; i<D/5; ++i){
        MatMul(Ar[5*i+1], Ai[5*i+1], Ar[5*i],     Ai[5*i],     Br[i], Bi[i]);
        MatMul(Ar[5*i+2], Ai[5*i+2], Br[i],       Bi[i],       Cr[i], Ci[i]);
        MatMul(Ar[5*i+3], Ai[5*i+3], Cr[i],       Ci[i],       Dr[i], Di[i]);
        MatMul(Ar[5*i+4], Ai[5*i+4], Dr[i],       Di[i],       Er[i], Ei[i]);
    }
    tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D/5, LOG_qdiv>(
        Er, Ei, ct_up, ct_down, skey, Act_up, Act_down
    );
}


template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int LOG_qdiv>
void tuple_grouped7_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    assert(D % 7 == 0);
    SparseDiagonal<1<<(LOGN-1), S* S>                   Br[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S>                   Bi[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S>                Cr[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S>                Ci[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S>             Dr[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S>             Di[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S>          Er[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S>          Ei[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S* S>       Fr[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S* S>       Fi[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S* S* S>    Gr[D/7];
    SparseDiagonal<1<<(LOGN-1), S* S* S* S* S* S* S>    Gi[D/7];

    for (int i=0; i<D/7; ++i){
        MatMul(Ar[7*i+1], Ai[7*i+1], Ar[7*i],     Ai[7*i],     Br[i], Bi[i]);
        MatMul(Ar[7*i+2], Ai[7*i+2], Br[i],       Bi[i],       Cr[i], Ci[i]);
        MatMul(Ar[7*i+3], Ai[7*i+3], Cr[i],       Ci[i],       Dr[i], Di[i]);
        MatMul(Ar[7*i+4], Ai[7*i+4], Dr[i],       Di[i],       Er[i], Ei[i]);
        MatMul(Ar[7*i+5], Ai[7*i+5], Er[i],       Ei[i],       Fr[i], Fi[i]);
        MatMul(Ar[7*i+6], Ai[7*i+6], Fr[i],       Fi[i],       Gr[i], Gi[i]);
    }
    tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D/7, LOG_qdiv>(
        Gr, Gi, ct_up, ct_down, skey, Act_up, Act_down
    );
}


template<int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G, int LOG_qdiv>
void tuple_grouped_serial_linear_transform_adaptive(
    const SparseDiagonal<1<<(LOGN-1), S> Ar[D],
    const SparseDiagonal<1<<(LOGN-1), S> Ai[D],
    const R_Q_square<LOGQ, 1<<LOGN>& ct_up,
    const R_Q_square<LOGQ, 1<<LOGN>& ct_down,
    const int skey[1<<LOGN],
    R_Q_square<LOGQ, 1<<LOGN>& Act_up,
    R_Q_square<LOGQ, 1<<LOGN>& Act_down
){
    if constexpr (G == 1){
        tuple_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(
            Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down
        );
    } else if constexpr (G == 2){
        tuple_grouped2_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(
            Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down
        );
    } else if constexpr (G == 3){
        tuple_grouped3_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(
            Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down
        );
    } else if constexpr (G == 4){
        tuple_grouped4_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(
            Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down
        );
    } else if constexpr (G == 5){
        tuple_grouped5_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(
            Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down
        );
    } else if constexpr (G == 7){
        tuple_grouped7_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D, LOG_qdiv>(
            Ar, Ai, ct_up, ct_down, skey, Act_up, Act_down
        );
    } else {
        static_assert(G==1 || G==2 || G==3 || G==4 || G==5 || G==7,
                      "Unsupported group size G");
    }
}

