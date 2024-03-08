// @file signaturecontext.cpp - Implementation file for signature context class
// @author TPOC: contact@palisade-crypto.org
//
// @copyright Copyright (c) 2019, New Jersey Institute of Technology (NJIT)
// All rights reserved.
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution. THIS SOFTWARE IS
// PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
// EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "signaturecontext.h"

namespace lbcrypto {
// Method for setting up a GPV context with specific parameters
template <class Element>
void SignatureContext<Element>::GenerateGPVContext(usint ringsize, usint bits,
                                                   usint base, bool VerifyNorm, usint dimension, usint k0, usint k1) {
  usint sm = ringsize * 2;
  double stddev = 0;
  if (k0 == 1 && k1 == 1){
      stddev = SIGMA;
  } else {
      stddev = ETA;
  }
  typename Element::DggType dgg(stddev);
  typename Element::Integer smodulus;
  typename Element::Integer srootOfUnity;

  smodulus = FirstPrime<typename Element::Integer>(bits, sm);
  srootOfUnity = RootOfUnity(sm, smodulus);
  ILParamsImpl<typename Element::Integer> ilParams =
      ILParamsImpl<typename Element::Integer>(sm, smodulus, srootOfUnity);

  ChineseRemainderTransformFTT<typename Element::Vector>::PreCompute(
      srootOfUnity, sm, smodulus);
  
  DiscreteFourierTransform::PreComputeTable(sm);

  auto silparams = std::make_shared<ILParamsImpl<typename Element::Integer>>(ilParams);
  m_params = std::make_shared<GPVSignatureParameters<Element>>(silparams, dgg, base, dimension, VerifyNorm, k0, k1);
  m_scheme = std::make_shared<GPVSignatureScheme<Element>>();
}
// Method for setting up a GPV context with desired security level only
template <class Element>
void SignatureContext<Element>::GenerateGPVContext(usint ringsize, bool VerifyNorm) {
  usint base, k, k0, k1, dimension;
  switch (ringsize) {
//    case 16:
//      k = 2;
//      base = 8;
    case 256:
      k = 52;
      k0 = 6;
      k1 = 13;
      base = 8;
      dimension = 4;
      GenerateGPVContext(ringsize, k, base, VerifyNorm, dimension, k0, k1);
      break;
    case 512:
      k = 13;
      k0 = 7;
      k1 = 15;
      base = 8;
      dimension = 4;
      GenerateGPVContext(ringsize, k, base, VerifyNorm, dimension, k0, k1);
      break;
    case 1024:
      k = 27;
      base = 64;
      GenerateGPVContext(ringsize, k, base, VerifyNorm);
      break;
    default:
      PALISADE_THROW(config_error, "Unknown ringsize");
  }

}
// Method for key generation
template <class Element>
void SignatureContext<Element>::KeyGen(LPSignKey<Element>* sk, LPVerificationKey<Element>* vk, bool setup) {
  m_scheme->KeyGen(m_params, sk, vk, setup);
}

// Method for signing a given plaintext
template <class Element>
void SignatureContext<Element>::Sign(const LPSignPlaintext<Element>& pt,
                                     const LPSignKey<Element>& sk,
                                     const LPVerificationKey<Element>& vk,
                                     LPSignature<Element>* sign) {
  m_scheme->Sign(m_params, sk, vk, pt, sign);
}

// Method for generate CRS
template <class Element>
void SignatureContext<Element>::CrsGen(const LPVerificationKey<Element>& vki,
                                       const LPSignKey<Element>& sk,
                                       const LPVerificationKey<Element>& vk,
                                       const LPSignPlaintext<Element>& seed,
                                       LPSignature<Element>* sign) {
    m_scheme->CrsGen(m_params, sk, vk, vki, seed, sign);
}

// Method for offline phase of signing a given plaintext
template <class Element>
void SignatureContext<Element>::SignOfflinePhase(
    const LPSignKey<Element>& sk, PerturbationVector<Element>& pv) {
  pv = m_scheme->SampleOffline(m_params, sk);
}
// Method for online phase of signing a given plaintext
template <class Element>
void SignatureContext<Element>::SignOnlinePhase(
    const LPSignPlaintext<Element>& pt, const LPSignKey<Element>& sk,
    const LPVerificationKey<Element>& vk, const PerturbationVector<Element> pv,
    LPSignature<Element>* signatureText) {
  m_scheme->SignOnline(m_params, sk, vk, pv, pt, signatureText);
}
// Method for verifying the plaintext and signature
template <class Element>
bool SignatureContext<Element>::Verify(const LPSignPlaintext<Element>& pt,
                                       const LPSignature<Element>& signature,
                                       const LPVerificationKey<Element>& vk,
                                       const LPSignPlaintext<Element>& seed) {
  return m_scheme->Verify(m_params, vk, signature, pt, seed);
}

// Method for verifying the plaintext and signature
template <class Element>
bool SignatureContext<Element>::VerifyMulti(const LPSignPlaintext<Element>& pt,
                                            const LPSignature<Element>& signature,
                                            const LPVerificationKey<Element>& vk,
                                            const Matrix<Element>& weight,
                                            const LPSignPlaintext<Element> seeds[]) {
    return m_scheme->VerifyMulti(m_params, vk, signature, pt, weight, seeds);
}

// Method for verifying the plaintext and signature
template <class Element>
bool SignatureContext<Element>::VerifyAggre(const LPSignPlaintext<Element> pts[],
                                            const LPSignature<Element>& signature,
                                            const LPVerificationKey<Element>& vk,
                                            const Matrix<Element>& weight,
                                            const LPSignPlaintext<Element> seeds[]) {
    return m_scheme->VerifyAggre(m_params, vk, signature, pts, weight, seeds);
}

}  // namespace lbcrypto
