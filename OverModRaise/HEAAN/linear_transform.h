#pragma once

#include "core.h"
#include "endecode.h"
#include <vector>    
#include <map>        
#include <algorithm>  
#include <set>       
#include <cmath>    

template<int N>
void rot( const int s[N], int s_rot[N] ){
	for(int i=0; i<N; i++){
		int q = (5*i)/N;
		int r = (5*i)%N;
		s_rot[r] = (q%2==0)?s[i]:-s[i];
	}
}

template<int N>
void rot( const int s[N], int s_rot[N], const int r){
	// pow = (5^r) % (2*N)
	int pow = 1;
	for(int i = 0; i < r; ++i) {
		pow = (5 * pow) % (2*N);
	}
	int j = 0;
	for(int i = 0; i < N; ++i, j = (j + pow) % (2*N)) {
		if(j < N)
			s_rot[j] = s[i];
		else
			s_rot[j-N] = -s[i];
	}
}

template<int LOGQ, int N>
void rot( const R_Q<LOGQ,N>& pt, 
				R_Q<LOGQ,N>& pt_rot, const int r = 1){
	int pow = 1;
	for(int i = 0; i < r; ++i) {
		pow = (5 * pow) % (2*N);
	}
	int j = 0;
	for(int i=0; i<N; i++, j = (j + pow) % (2*N)){
		if(j < N)
			pt_rot[j] = pt[i];
		else {
			pt_rot[j-N] = pt[i];
			pt_rot[j-N].negate();
		}
	}
}

//-------------------------------------------------------------------------------
// rotate ct off times
//-------------------------------------------------------------------------------
template<int LOGQ, int N>
void rot_ct(const R_Q_square<  LOGQ, N>& ct, int off,
	        const R_Q_square<2*LOGQ, N>& rkey,
	              R_Q_square<  LOGQ, N>& ct_rot) {
	rot<LOGQ, N>(ct[0], ct_rot[0], off);
	rot<LOGQ, N>(ct[1], ct_rot[1], off);
	HEAAN<LOGQ, N>::ks(rkey, ct_rot);
}

template<int N>
void conj( const int s[N], int s_conj[N] ){
	s_conj[0] = s[0];
	for(int i=1; i<N; i++) s_conj[i] = -s[N-i];
}

template<int LOGQ, int N>
void conj( const R_Q<LOGQ,N>& pt,
				 R_Q<LOGQ,N>& pt_conj ){
	pt_conj[0] = pt[0];
	for(int i=1; i<N; i++){
		pt_conj[i] = pt[N-i];
		pt_conj[i].negate();
	}
}

template<int LOGQ, int N>
void conj( const R_Q_square<  LOGQ,N>& ct,
		   const R_Q_square<2*LOGQ,N>& ckey,
				 R_Q_square<  LOGQ,N>& ct_conj ){
	conj<LOGQ,N>(ct[0], ct_conj[0]);
	conj<LOGQ,N>(ct[1], ct_conj[1]);
	HEAAN<LOGQ,N>::ks(ckey,ct_conj);
}

template<int N>
void mat_to_vecs( const double A[N][N],
						double V[N][N]){
	for(int i=0; i<N; i++)
	for(int j=0; j<N; j++)
		V[i][j] = A[j][(i+j)%N];
}

// naive. not considering sparse. one slot rot by one for loop
template<int LOGQ, int LOGN, int LOGDELTA>
void linear_transform( const double Ar[1 << (LOGN - 1)][1 << (LOGN - 1)],
					   const double Ai[1 << (LOGN - 1)][1 << (LOGN - 1)],
					   const R_Q_square<  LOGQ, 1 << LOGN>& ct,
					   const R_Q_square<2*LOGQ, 1 << LOGN>& rkey,
							 R_Q_square<  LOGQ, 1 << LOGN>& Act ){
	const int N = 1 << LOGN;

	double vr[N/2][N/2]; mat_to_vecs<N/2>(Ar,vr);
	double vi[N/2][N/2]; mat_to_vecs<N/2>(Ai,vi);
	Act.setzero(); R_Q_square<LOGQ,N> ct_rot = ct;
	for(int i=0; i<N/2; i++){
		R_Q<LOGQ,N> pt;
		encode<LOGQ, LOGN>(vr[i],vi[i], 1ULL<<LOGDELTA, pt);		

		R_Q_square<LOGQ,N> ct_rot_pt = ct_rot;
		ct_rot_pt *= pt;

		Act += ct_rot_pt; R_Q_square<LOGQ,N> temp; 
		rot_ct<LOGQ,N>(ct_rot, rkey, temp); // BUG
		ct_rot = temp; 
	}
}


//-----------------------------------------------------------------------------------
// A : sparse diagonal matrix
// Act: lineartransform(A,ct) such that
// decode(dec(Act,s,Q),Delta^2) = A. decode(dec(ct,s,Q),Delta)
//-----------------------------------------------------------------------------------
template< int LOGQ, int LOGN, int LOGDELTA, int S >
void linear_transform( const SparseDiagonal<1 << (LOGN - 1),S>& Ar,
					   const SparseDiagonal<1 << (LOGN - 1),S>& Ai,
					   const R_Q_square<  LOGQ, 1 << LOGN>& ct,
					   const R_Q_square<2*LOGQ, 1 << LOGN>* rkey[S],
							 R_Q_square<  LOGQ, 1 << LOGN>& Act ){
	const int N = 1 << LOGN;
	Act.setzero();
	for (int s = 0; s < S; s++) {
		if (Ar.zero[s] == false || Ai.zero[s] == false) {
			R_Q<LOGQ, N> pt;
			encode<LOGQ, LOGN>(Ar.vec[s],
							Ai.vec[s], 1ULL << LOGDELTA, pt);
			R_Q_square<LOGQ, N> ct_rot;
			if (Ar.off[s] != 0)
				rot_ct<LOGQ, N>(ct, Ar.off[s],*(rkey[s]), ct_rot);
			ct_rot *= pt;
			Act += ct_rot;
		}
	}
}


template< int LOGQ, int LOGN, int LOGDELTA, int S >
void linear_transform(const SparseDiagonal<1 << (LOGN - 1), S>& Ar,
	const SparseDiagonal<1 << (LOGN - 1), S>& Ai,
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
	const int N = 1 << LOGN;
	Act.setzero();
	for (int s = 0; s < S; s++) {
		assert(Ar.off[s] == Ai.off[s]);
		if (Ar.zero[s] && Ai.zero[s])
			continue;
		R_Q<LOGQ, N> pt;
		encode<LOGQ, LOGN>(Ar.vec[s],
			Ai.vec[s], 1ULL << LOGDELTA, pt);
		R_Q_square<LOGQ, N> ct_rot(ct);
		if (Ar.off[s] != 0) {
			R_Q_square<2 * LOGQ, 1 << LOGN> rkey;
			int skey_rot[N];
			rot<N>(skey, skey_rot, Ar.off[s]);
			HEAAN<LOGQ, N>::swkgen(skey_rot, skey, rkey);
			rot_ct<LOGQ, N>(ct, Ar.off[s], rkey, ct_rot);
		}
		ct_rot *= pt;
		Act += ct_rot;
	}
}

















template<int LOGQ, int LOGN, int LOGDELTA, int S>
void linear_transform_sparse_babystep_giantstep(
    const SparseDiagonal<(1<<(LOGN-1)), S>& Ar,
    const SparseDiagonal<(1<<(LOGN-1)), S>& Ai,
    const R_Q_square<LOGQ, (1<<LOGN)>& ct,
    const int skey[(1<<LOGN)],
    R_Q_square<LOGQ, (1<<LOGN)>& Act,
    int baby_step_size 
){
    constexpr int N    = 1 << LOGN;
    constexpr int HALF = N >> 1;

    auto normH = [](int x)->int{ int y = x % HALF; return (y < 0) ? (y + HALF) : y; };
    auto igcd  = [](int a, int b){ a = std::abs(a); b = std::abs(b);
                                   while (b){ int t = a % b; a = b; b = t; }
                                   return a; };

    auto get_off = [&](int s)->int {
        if (!Ar.zero[s] && !Ai.zero[s]) {
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
    Act.setzero();
    for (int s=0; s<S; ++s) {
        if (!Ar.zero[s] || !Ai.zero[s]) {
            offs.push_back(normH(get_off(s)));
        }
    }
    if (offs.empty()) return;

    int k = 0;
    for (int v : offs) if (v) k = (k==0 ? v : igcd(k, v));
    if (k == 0) k = 1;
    std::cerr << "[BSGS] stride k = " << k << '\n';  

    for (int v : offs) { assert(normH(v) % k == 0 && "offset not multiple of stride k"); }

  
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

   
    std::map<int, R_Q_square<LOGQ, N>> ct_baby; 
    for (int i : babyI){
        int amt = normH(i * k);
        if (amt == 0) { ct_baby[i] = ct; continue; }
        R_Q_square<2*LOGQ, N> rkey_i;
        int skey_rot_i[N];
        rot<N>(skey, skey_rot_i, amt);
        HEAAN<LOGQ, N>::swkgen(skey_rot_i, skey, rkey_i);
        rot_ct<LOGQ, N>(ct, amt, rkey_i, ct_baby[i]);
    }
    if (!ct_baby.count(0)) ct_baby[0] = ct;

  
    for (auto &kv : groupJ){
        const int j = kv.first;
        const auto& idxs = kv.second;

        R_Q_square<LOGQ, N> inner; inner.setzero();

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

            R_Q<LOGQ, N> pt;
            encode<LOGQ, LOGN>(re.data(), im.data(), 1ULL<<LOGDELTA, pt);

            
            R_Q_square<LOGQ, N> tmp = ct_baby[i];
            tmp *= pt;
            inner += tmp;
        }

       
        if (j != 0){
            int amt = r; 
            R_Q_square<2*LOGQ, N> rkey_j;
            int skey_rot_j[N];
            rot<N>(skey, skey_rot_j, amt);
            HEAAN<LOGQ, N>::swkgen(skey_rot_j, skey, rkey_j);

            R_Q_square<LOGQ, N> inner_rot;
            rot_ct<LOGQ, N>(inner, amt, rkey_j, inner_rot);
            Act += inner_rot;
        } else {
            Act += inner;
        }
    }
}


// ============================================================================
// Adaptive wrapper 
// ============================================================================
template< int LOGQ, int LOGN, int LOGDELTA, int S >
void linear_transform_adaptive_sparse_babystep_giantstep(
    const SparseDiagonal<1 << (LOGN - 1), S>& Ar,
    const SparseDiagonal<1 << (LOGN - 1), S>& Ai,
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act
){
    int cnt = 0;
    for (int s = 0; s < S; ++s) if (!Ar.zero[s] || !Ai.zero[s]) ++cnt;
    int suggested_baby = std::max(1, (int)std::sqrt((double)cnt));

    linear_transform_sparse_babystep_giantstep<LOGQ, LOGN, LOGDELTA, S>(
        Ar, Ai, ct, skey, Act, suggested_baby
    );
}

























template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
	R_Q_square<LOGQ, 1 << LOGN> ct_temp;
	Act = ct;
	for (int d = 0; d < D; ++d) {
		ct_temp = Act;
		linear_transform<LOGQ, LOGN, LOGDELTA, S>(Ar[d], Ai[d], ct_temp, skey, Act);
	}
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped2_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN >& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
	if (D % 2 == 0) {
		SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
		SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

		for (int i = 0; i < D / 2; ++i) {
			MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
		}
		serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, Act);
	}

	else {
		SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
		SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

		for (int i = 0; i < D / 2; ++i) {
			MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
		}
		R_Q_square<LOGQ, 1 << LOGN> ct_temp;
		serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, ct_temp);
		linear_transform<LOGQ, LOGN, LOGDELTA, S>(Ar[D-1], Ai[D-1], ct_temp, skey, Act);
	}

	
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped3_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
	assert(D % 3 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 3];

	for (int i = 0; i < D / 3; ++i) {
		MatMul(Ar[3 * i + 1], Ai[3 * i + 1], Ar[3 * i], Ai[3 * i], Br[i], Bi[i]);
		MatMul(Ar[3 * i + 2], Ai[3 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
	}
	serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S, D / 3>(Cr, Ci, ct, skey, Act);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped4_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
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
	serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4>(Dr, Di, ct, skey, Act);
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

	R_Q_square<LOGQ, 1 << LOGN> ct_temp;
	serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4>(Dr, Di, ct, skey, ct_temp);
	linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S>(LCr, LCi, ct_temp, skey, Act);
}
}


template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped5_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
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
	serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D / 5>(Er, Ei, ct, skey, Act);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped7_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
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
	serial_linear_transform<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D / 7>(Gr, Gi, ct, skey, Act);
}




template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G = 1>
void grouped_serial_linear_transform(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
	if (G == 1)
		serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 2)
		grouped2_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 3)
		grouped3_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 4)
		grouped4_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 5)
		grouped5_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 7)
		grouped7_serial_linear_transform<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else
		assert(false);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S >
void linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S>& Ar,
	const SparseDiagonal<1 << (LOGN - 1), S>& Ai,
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
	const int N = 1 << LOGN;
	Act.setzero();
	for (int s = 0; s < S; s++) {
		assert(Ar.off[s] == Ai.off[s]);
		if (Ar.zero[s] && Ai.zero[s])
			continue;
		R_Q<LOGQ, N> pt;
		R_Q<LOGQ, N> pt_rot;
		encode<LOGQ, LOGN>(Ar.vec[s], Ai.vec[s], 1ULL << LOGDELTA, pt);

		R_Q_square<LOGQ, N> ct_rot(ct);
		R_Q_square<LOGQ - LOGDELTA, N> ct_block;

		if (Ar.off[s] != 0) {

			rot<LOGQ, N>(pt, pt_rot, N/2 - Ar.off[s]);
			ct_rot = ct;
			ct_rot *= pt_rot;
			R_Q_square<LOGQ - LOGDELTA, N> ct_rs;
			

			RS<LOGQ, LOGQ - LOGDELTA, N>(ct_rot, ct_rs);

			R_Q_square<2 * (LOGQ - LOGDELTA), 1 << LOGN> rkey;
			int skey_rot[N];
			rot<N>(skey, skey_rot, Ar.off[s]);
			HEAAN<LOGQ - LOGDELTA, N>::swkgen(skey_rot, skey, rkey);
			rot_ct<LOGQ - LOGDELTA, N>(ct_rs, Ar.off[s], rkey, ct_block);
		}

		if (Ar.off[s] == 0) {
			ct_rot = ct;
			ct_rot *= pt;
			RS<LOGQ, LOGQ - LOGDELTA, N>(ct_rot, ct_block);
		}
		
		Act += ct_block;
	}
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
	
	Act.setzero();

	linear_transform_sw<LOGQ, LOGN, LOGDELTA, S>(Ar[0], Ai[0], ct, skey, Act);

	R_Q_square<LOGQ - LOGDELTA, 1 << LOGN> ct_temp;

	for (int d = 1; d < D; ++d) {
		ct_temp = Act;
		linear_transform<LOGQ - LOGDELTA, LOGN, LOGDELTA, S>(Ar[d], Ai[d], ct_temp, skey, Act);
	}
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped2_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ -  LOGDELTA, 1 << LOGN>& Act) {
	if(D % 2 == 0){
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

	for (int i = 0; i < D / 2; ++i) {
		MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
	}
	serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, Act);
}

else {
		SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
		SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

		for (int i = 0; i < D / 2; ++i) {
			MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
		}
		R_Q_square<LOGQ -  LOGDELTA, 1 << LOGN> ct_temp;
		serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, ct_temp);
		linear_transform<LOGQ -  LOGDELTA, LOGN, LOGDELTA, S>(Ar[D-1], Ai[D-1], ct_temp, skey, Act);
	}

}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped3_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
	assert(D % 3 == 0);
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 3];
	SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 3];

	for (int i = 0; i < D / 3; ++i) {
		MatMul(Ar[3 * i + 1], Ai[3 * i + 1], Ar[3 * i], Ai[3 * i], Br[i], Bi[i]);
		MatMul(Ar[3 * i + 2], Ai[3 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
	}
	serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S, D / 3>(Cr, Ci, ct, skey, Act);
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped4_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ -  LOGDELTA, 1 << LOGN>& Act) {
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
	serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4>(Dr, Di, ct, skey, Act);
}

	if(D % 4== 3){
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

	R_Q_square<LOGQ - LOGDELTA, 1 << LOGN> ct_temp;
	serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4>(Dr, Di, ct, skey, ct_temp);
	linear_transform<LOGQ - LOGDELTA, LOGN, LOGDELTA, S* S* S>(LCr, LCi, ct_temp, skey, Act);
}

}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped5_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
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
	serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D / 5>(Er, Ei, ct, skey, Act);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped7_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
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
	serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D / 7>(Gr, Gi, ct, skey, Act);
}


template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G > 
void grouped_serial_linear_transform_sw(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
	if (G == 1)
		serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 2)
		grouped2_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 3)
		grouped3_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 4)
		grouped4_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 5)
		grouped5_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else if (G == 7)
		grouped7_serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
	else
		assert(false);
}








template< int LOGQ, int LOGN, int LOGDELTA, int S, int D >
void serial_linear_transform_adaptive(
    const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act);

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G >
void grouped_serial_linear_transform_adaptive(
    const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act);







template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act) {
    
    R_Q_square<LOGQ, 1 << LOGN> ct_temp;
    Act = ct;
    for (int d = 0; d < D; ++d) {
        ct_temp = Act;
        linear_transform_adaptive_sparse_babystep_giantstep<LOGQ, LOGN, LOGDELTA, S>(
            Ar[d], Ai[d], ct_temp, skey, Act
        );
    }
}




template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped4_serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act) {

    if (D == 0) { Act.setzero(); return; }

    if (D % 4 == 0) {
        SparseDiagonal<1 << (LOGN - 1), S* S>       Br[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S* S>       Bi[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S* S* S>    Cr[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S* S* S>    Ci[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[D / 4];

        for (int i = 0; i < D / 4; ++i) {
            MatMul(Ar[4 * i + 1], Ai[4 * i + 1], Ar[4 * i],     Ai[4 * i],     Br[i], Bi[i]);
            MatMul(Ar[4 * i + 2], Ai[4 * i + 2], Br[i],         Bi[i],         Cr[i], Ci[i]);
            MatMul(Ar[4 * i + 3], Ai[4 * i + 3], Cr[i],         Ci[i],         Dr[i], Di[i]);
        }
        serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S, D / 4>(
            Dr, Di, ct, skey, Act);
        return;
    }

    if (D % 4 == 3) {
        SparseDiagonal<1 << (LOGN - 1), S* S>    LBr;
        SparseDiagonal<1 << (LOGN - 1), S* S>    LBi;
        SparseDiagonal<1 << (LOGN - 1), S* S* S> LCr;
        SparseDiagonal<1 << (LOGN - 1), S* S* S> LCi;

        MatMul(Ar[1], Ai[1], Ar[0], Ai[0], LBr, LBi);
        MatMul(Ar[2], Ai[2], LBr,  LBi,  LCr, LCi);

        R_Q_square<LOGQ, 1 << LOGN> ct_after3;
        linear_transform_adaptive_sparse_babystep_giantstep<LOGQ, LOGN, LOGDELTA, S* S* S>(
            LCr, LCi, ct, skey, ct_after3);

        if (D == 3) { Act = ct_after3; return; } 

        constexpr int R = (D - 3) / 4; 
        SparseDiagonal<1 << (LOGN - 1), S* S>       Br[R];
        SparseDiagonal<1 << (LOGN - 1), S* S>       Bi[R];
        SparseDiagonal<1 << (LOGN - 1), S* S* S>    Cr[R];
        SparseDiagonal<1 << (LOGN - 1), S* S* S>    Ci[R];
        SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Dr[R];
        SparseDiagonal<1 << (LOGN - 1), S* S* S* S> Di[R];

        for (int i = 0; i < R; ++i) {
            const int base = 3 + 4 * i;
            MatMul(Ar[base + 1], Ai[base + 1], Ar[base + 0], Ai[base + 0], Br[i], Bi[i]);
            MatMul(Ar[base + 2], Ai[base + 2], Br[i],        Bi[i],        Cr[i], Ci[i]);
            MatMul(Ar[base + 3], Ai[base + 3], Cr[i],        Ci[i],        Dr[i], Di[i]);
        }

        serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S, R>(
            Dr, Di, ct_after3, skey, Act);
        return;
    }

    assert(false && "grouped4_serial_linear_transform_adaptive: D % 4 must be 0 or 3");
}




// G = 2
template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped2_serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act) {

    assert(D % 2 == 0);

    SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
    SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

    for (int i = 0; i < D / 2; ++i) {
        MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
    }
    serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, Act);
}

// G = 3
template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped3_serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act) {

    assert(D % 3 == 0);

    SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 3];
    SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 3];
    SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 3];
    SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 3];

    for (int i = 0; i < D / 3; ++i) {
        MatMul(Ar[3 * i + 1], Ai[3 * i + 1], Ar[3 * i], Ai[3 * i], Br[i], Bi[i]);
        MatMul(Ar[3 * i + 2], Ai[3 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
    }
    serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S, D / 3>(Cr, Ci, ct, skey, Act);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped5_serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ, 1 << LOGN>& Act) {
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
	serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D / 5>(Er, Ei, ct, skey, Act);
}

// G = 7
template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped7_serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act) {

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
    serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D / 7>(Gr, Gi, ct, skey, Act);
}


template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G>
void grouped_serial_linear_transform_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ, 1 << LOGN>& Act) {

    if (G == 1)
        serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 2)
        grouped2_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 3)
        grouped3_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 4)
        grouped4_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 5)
        grouped5_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 7)
        grouped7_serial_linear_transform_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else
        assert(false);
}







template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
	
	Act.setzero();

	linear_transform_sw<LOGQ, LOGN, LOGDELTA, S>(Ar[0], Ai[0], ct, skey, Act);

	R_Q_square<LOGQ - LOGDELTA, 1 << LOGN> ct_temp;

	for (int d = 1; d < D; ++d) {
		ct_temp = Act;
		linear_transform_adaptive_sparse_babystep_giantstep<LOGQ - LOGDELTA, LOGN, LOGDELTA, S>(Ar[d], Ai[d], ct_temp, skey, Act);
	}
}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped2_serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ -  LOGDELTA, 1 << LOGN>& Act) {
	if(D % 2 == 0){
	SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
	SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

	for (int i = 0; i < D / 2; ++i) {
		MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
	}
	serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, Act);
}

else {
		SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 2];
		SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 2];

		for (int i = 0; i < D / 2; ++i) {
			MatMul(Ar[2 * i + 1], Ai[2 * i + 1], Ar[2 * i], Ai[2 * i], Br[i], Bi[i]);
		}
		R_Q_square<LOGQ -  LOGDELTA, 1 << LOGN> ct_temp;
		serial_linear_transform_sw<LOGQ, LOGN, LOGDELTA, S* S, D / 2>(Br, Bi, ct, skey, ct_temp);
		linear_transform<LOGQ -  LOGDELTA, LOGN, LOGDELTA, S>(Ar[D-1], Ai[D-1], ct_temp, skey, Act);
	}

}



template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped4_serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ - LOGDELTA, 1 << LOGN>& Act) 
{
    if (D == 0) { Act.setzero(); return; }

    if (D % 4 == 0) {
        SparseDiagonal<1 << (LOGN - 1), S*S> Br[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S*S> Bi[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S*S*S> Cr[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S*S*S> Ci[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S*S*S*S> Dr[D / 4];
        SparseDiagonal<1 << (LOGN - 1), S*S*S*S> Di[D / 4];

        for (int i = 0; i < D / 4; ++i) {
            MatMul(Ar[4*i + 1], Ai[4*i + 1], Ar[4*i], Ai[4*i], Br[i], Bi[i]);
            MatMul(Ar[4*i + 2], Ai[4*i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
            MatMul(Ar[4*i + 3], Ai[4*i + 3], Cr[i], Ci[i], Dr[i], Di[i]);
        }
        serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S*S*S*S, D / 4>(Dr, Di, ct, skey, Act);
    }
    else if (D % 4 == 3) {
        if (D == 3) {
            SparseDiagonal<1 << (LOGN - 1), S*S> Br;
            SparseDiagonal<1 << (LOGN - 1), S*S> Bi;
            SparseDiagonal<1 << (LOGN - 1), S*S*S> Cr;
            SparseDiagonal<1 << (LOGN - 1), S*S*S> Ci;
            MatMul(Ar[1], Ai[1], Ar[0], Ai[0], Br, Bi);
            MatMul(Ar[2], Ai[2], Br, Bi, Cr, Ci);
            linear_transform_sw<LOGQ, LOGN, LOGDELTA, S*S*S>(Cr, Ci, ct, skey, Act);
        } else {
            R_Q_square<LOGQ - LOGDELTA, 1 << LOGN> ct_temp;
            SparseDiagonal<1 << (LOGN - 1), S*S> First_Br;
            SparseDiagonal<1 << (LOGN - 1), S*S> First_Bi;
            SparseDiagonal<1 << (LOGN - 1), S*S*S> First_Cr;
            SparseDiagonal<1 << (LOGN - 1), S*S*S> First_Ci;
            MatMul(Ar[1], Ai[1], Ar[0], Ai[0], First_Br, First_Bi);
            MatMul(Ar[2], Ai[2], First_Br, First_Bi, First_Cr, First_Ci);
            linear_transform_sw<LOGQ, LOGN, LOGDELTA, S*S*S>(First_Cr, First_Ci, ct, skey, ct_temp);

            grouped_serial_linear_transform_adaptive<LOGQ - LOGDELTA, LOGN, LOGDELTA, S, D - 3, 4>(
                &Ar[3], &Ai[3], ct_temp, skey, Act
            );
        }
    }
}


// -----------------------------------------------------------------------------
//  G = 3
// -----------------------------------------------------------------------------
template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped3_serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ - LOGDELTA, 1 << LOGN>& Act) 
{
    static_assert(D > 0, "D must be positive");
    assert(D % 3 == 0);

    SparseDiagonal<1 << (LOGN - 1), S* S> Br[D / 3];
    SparseDiagonal<1 << (LOGN - 1), S* S> Bi[D / 3];
    SparseDiagonal<1 << (LOGN - 1), S* S* S> Cr[D / 3];
    SparseDiagonal<1 << (LOGN - 1), S* S* S> Ci[D / 3];

    for (int i = 0; i < D / 3; ++i) {
        MatMul(Ar[3 * i + 1], Ai[3 * i + 1], Ar[3 * i], Ai[3 * i], Br[i], Bi[i]);
        MatMul(Ar[3 * i + 2], Ai[3 * i + 2], Br[i], Bi[i], Cr[i], Ci[i]);
    }
    serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S, D / 3>(Cr, Ci, ct, skey, Act);
}

template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped5_serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
	const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
	const R_Q_square<  LOGQ, 1 << LOGN>& ct,
	const int skey[1 << LOGN],
	R_Q_square<  LOGQ - LOGDELTA, 1 << LOGN>& Act) {
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
	serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S* S, D / 5>(Er, Ei, ct, skey, Act);
}

// -----------------------------------------------------------------------------
// G = 7
// -----------------------------------------------------------------------------
template< int LOGQ, int LOGN, int LOGDELTA, int S, int D>
void grouped7_serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ - LOGDELTA, 1 << LOGN>& Act) 
{
    static_assert(D > 0, "D must be positive");
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
    serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S* S* S* S* S* S* S, D / 7>(Gr, Gi, ct, skey, Act);
}

// -----------------------------------------------------------------------------
// Dispatcher (G = 1,2,3,4,5,7) for _sw_adaptive (skey only)
// -----------------------------------------------------------------------------
template< int LOGQ, int LOGN, int LOGDELTA, int S, int D, int G >
void grouped_serial_linear_transform_sw_adaptive(const SparseDiagonal<1 << (LOGN - 1), S> Ar[D],
    const SparseDiagonal<1 << (LOGN - 1), S> Ai[D],
    const R_Q_square<LOGQ, 1 << LOGN>& ct,
    const int skey[1 << LOGN],
    R_Q_square<LOGQ - LOGDELTA, 1 << LOGN>& Act) 
{
    if (G == 1)
        serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 2)
        grouped2_serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 3)
        grouped3_serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 4)
        grouped4_serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 5)
        grouped5_serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else if (G == 7)
        grouped7_serial_linear_transform_sw_adaptive<LOGQ, LOGN, LOGDELTA, S, D>(Ar, Ai, ct, skey, Act);
    else
        assert(false);
}
