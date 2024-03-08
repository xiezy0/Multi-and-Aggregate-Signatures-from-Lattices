// @file
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
//
// @section DESCRIPTION
// This code provides the utility for GPV Ring-LWE signature scheme with
// trapdoors. The scheme implemented can be found in the paper
// https://eprint.iacr.org/2013/297.pdf. Construction 1 of the section 3.2 is
// used in this implementation.

#ifndef SIGNATURE_LWESIGN_H
#define SIGNATURE_LWESIGN_H

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "palisade.h"

 #include "signaturecore.h"

#include "lattice/trapdoor.h"

namespace lbcrypto {

/*
 *  @brief Templated class for holding signatures.
 *  @tparam is the ring element
 */
template <class Element>
class GPVSignature : public LPSignature<Element> {
 public:
  /**
   * Default constructor
   */
  GPVSignature() {}

  /**
   *Method for setting the element in signature
   *
   *@param signature Element vector to be set as the signature
   */
  void SetSignature(shared_ptr<Matrix<Element>> signature) {
    m_signature = signature;
  }
  /**
   *Method for getting the element in signature
   *
   *@return the element held as signature
   */
  const Matrix<Element>& GetSignature() const { return *m_signature; }

  /**
   *Destructor
   */
  ~GPVSignature() {}

 private:
  // Matrix of polynomials acting as actual signature
  shared_ptr<Matrix<Element>> m_signature;
  /*
   *@brief Overloaded dummy method
   */
  void forceImplement() {}
};

/*
 *  @brief Templated class for holding signatures.
 *  @tparam is the ring element
 */
template <class Element>
class GPVPlaintext : public LPSignPlaintext<Element> {
 public:
  /**
   * Default constructor
   */
  GPVPlaintext() {}
  /**
   *@brief Constructor
   *@param plaintext Plaintext to be set for signing
   */
  explicit GPVPlaintext(const string& plaintext) {
    this->m_plaintext = plaintext;
  }

  /**
   *Method for setting the element in plaintext
   *
   *@param plaintext Plaintext to be set for signing
   */
  void SetPlaintext(const string& plaintext) { this->m_plaintext = plaintext; }
  /**
   *Method for getting the element in plaintext
   *
   *@return the element held as plaintext
   */
  const string& GetPlaintext() const { return m_plaintext; }

  /**
   *Destructor
   */
  ~GPVPlaintext() {}
  string m_plaintext;
 private:
  // String as plaintext to be signed

  /*
   *@brief Overloaded dummy method
   */
  void forceImplement() {}
};

/**
 * @brief  Class holding parameters required for calculations in signature
 * schemes
 */
template <class Element>
class GPVSignatureParameters : public LPSignatureParameters<Element> {
 public:
  /**
   *Method for setting the ILParams held in this class
   *
   *@param params Parameters to be held, used in Element construction
   */
  void SetElemParams(shared_ptr<typename Element::Params> params,
                     usint base = 2, bool VerifyNormFlag = false) {
    m_params = params;
    m_base = base;
    const typename Element::Integer& q = params->GetModulus();
    size_t n = params->GetRingDimension();
    usint nBits = floor(log2(q.ConvertToDouble() - 1.0) + 1.0);
    m_k = ceil(nBits / log2(base));
    double c = (base + 1) * SIGMA;
    double s = SPECTRAL_BOUND(n, m_k, base);
    if (sqrt(s * s - c * c) <= KARNEY_THRESHOLD)
      m_dggLargeSigma = typename Element::DggType(sqrt(s * s - c * c));
    else
      m_dggLargeSigma = m_dgg;
    
    VerifyNorm = VerifyNormFlag;
  }

  /**
   *Method for accessing the ILParams held in this class
   *
   *@return Parameters held
   */
  shared_ptr<typename Element::Params> GetILParams() const { return m_params; }

  /**
   *Method for accessing the DiscreteGaussianGenerator object held in this class
   *
   *@return DiscreteGaussianGenerator object held
   */
  typename Element::DggType& GetDiscreteGaussianGenerator() { return m_dgg; }

  /**
   *Method for accessing the base for Gadget matrix
   *
   *@return the value of base held by the object
   */
  usint& GetBase() { return m_base; }


  usint& GetDimension() { return m_dimension; }

  usint& GetK0() { return m_k0; }
  usint& GetB0() { return m_b0; }

  usint& GetK1() { return m_k1; }
  usint& GetB1() { return m_b1; }

/**
   *Method for accessing the base for Gadget matrix
   *
   *@return the value of base held by the object
   */
  bool& GetVerifyNormFlag() { return VerifyNorm; }
  /**
   *Method for accessing the dimension for Gadget matrix
   *
   *@return the value of the dimension held by the object
   */
  usint& GetK() { return m_k; }

  /**
   *Method for accessing the DiscreteGaussianGenerator object held in this class
   *
   *@return DiscreteGaussianGenerator object held
   */
  typename Element::DggType& GetDiscreteGaussianGeneratorLargeSigma() {
    return m_dggLargeSigma;
  }

  typename Element::DggType& GetDiscreteGaussianGeneratorLargeSigma0() {
      return m_dggLargeSigma0;
  }

  typename Element::DggType& GetDiscreteGaussianGeneratorLargeSigma1() {
      return m_dggLargeSigma1;
  }

  /**
   *Constructor
   *@param params Parameters used in Element construction
   *@param dgg DiscreteGaussianGenerator used in sampling
   */
  GPVSignatureParameters(shared_ptr<typename Element::Params> params,
                         typename Element::DggType& dgg, usint base = 2, usint dimension = 1, bool VerifyNormFlag = false, usint k0 = 1, usint k1 = 1)
      : m_dgg(dgg), m_base(base), m_dimension(dimension), m_k0(k0), m_k1(k1) {
      const typename Element::Integer& q = params->GetModulus();
      size_t n = params->GetRingDimension();
      m_params = params;
      if (k0 == 1 && k1 == 1) {
          usint nBits = floor(log2(q.ConvertToDouble() - 1.0) + 1.0);
          //std::cout << nBits << std::endl;
          m_k = ceil(nBits / log2(base));
          double c = (base + 1) * SIGMA;
          double s = SPECTRAL_BOUND(n, m_k, base);
          if (sqrt(s * s - c * c) <= KARNEY_THRESHOLD)
              m_dggLargeSigma = typename Element::DggType(sqrt(s * s - c * c));
          else
              m_dggLargeSigma = m_dgg;
      } else {
          // b0 = k0 √q
          m_b0 = std::ceil(std::pow(q.ConvertToDouble(), 1.0/k0));
          m_b1 = std::ceil(std::pow(q.ConvertToDouble(), 1.0/k1));
//          std::cout << m_b0 << std::endl;
//          std::cout << m_b1 << std::endl;
          double S_0 = (1 + std::sqrt(m_k0/2)) * std::sqrt(dimension * n);
          double S_1 = (1 + std::sqrt(m_k1/2)) * std::sqrt(dimension * n);
          double sigmat0 = (S_0 + 1) * std::sqrt(std::pow(m_b0, 2) + 1) * ETA;
          double sigmat1 = (S_1 + 1) * std::sqrt(std::pow(m_b1, 2) + 1) * ETA;
          double s0 = SPECTRAL_BOUND_D_MP(n, m_k0, sigmat0, dimension);
          double s1 = SPECTRAL_BOUND_D_MP(n, m_k1, sigmat1, dimension);
          if (sqrt(s0 * s0 - sigmat0 * sigmat0) <= KARNEY_THRESHOLD)
              m_dggLargeSigma0 = typename Element::DggType(sqrt(s0 * s0 - sigmat0 * sigmat0));
          else
              m_dggLargeSigma0 = m_dgg;
          if (sqrt(s1 * s1 - sigmat1 * sigmat1) <= KARNEY_THRESHOLD)
              m_dggLargeSigma1 = typename Element::DggType(sqrt(s1 * s1 - sigmat1 * sigmat1));
          else
              m_dggLargeSigma1 = m_dgg;
      }
      VerifyNorm = VerifyNormFlag;
  }


//  GPVSignatureParameters(shared_ptr<typename Element::Params> params,
//                         typename Element::DggType& dgg, usint base = 2, usint dimension = 1, bool VerifyNormFlag = false, usint k0 = 1, usint k1 = 1)
//            : m_dgg(dgg), m_base(base), m_dimension(dimension), m_k0(k0), m_k1(k1) {
//        m_params = params;
//        const typename Element::Integer& q = params->GetModulus();
//        size_t n = params->GetRingDimension();
//
//        // b0 = k0 √q
//        m_b0 = std::ceil(std::pow(q.ConvertToDouble(), 1.0/k0));
//        m_b1 = std::ceil(std::pow(q.ConvertToDouble(), 1.0/k1));
//
//        double c0 = (m_b0 + 1) * MPSIGMA;
//        double c1 = (m_b1 + 1) * MPSIGMA;
//        double s0 = SPECTRAL_BOUND_D(n, m_k0, m_b0, dimension);
//        double s1 = SPECTRAL_BOUND_D(n, m_k1, m_b1, dimension);
//        if (sqrt(s0 * s0 - c0 * c0) <= KARNEY_THRESHOLD)
//            m_dggLargeSigma0 = typename Element::DggType(sqrt(s0 * s0 - c0 * c0));
//        else
//            m_dggLargeSigma0 = m_dgg;
//
//      if (sqrt(s1 * s1 - c1 * c1) <= KARNEY_THRESHOLD)
//          m_dggLargeSigma1 = typename Element::DggType(sqrt(s1 * s1 - c1 * c1));
//      else
//          m_dggLargeSigma1 = m_dgg;
//        VerifyNorm = VerifyNormFlag;
//
//        //m_b0 =
//    }

 private:
  // Parameters related to elements
  shared_ptr<typename Element::Params> m_params;
  // Discrete Gaussian Generator for random number generation
  typename Element::DggType m_dgg;
  // Discrete Gaussian Generator with high distribution parameter for random
  // number generation
  typename Element::DggType m_dggLargeSigma, m_dggLargeSigma0, m_dggLargeSigma1;
  // Trapdoor base
  usint m_base;
  // Trapdoor length
  usint m_k{};
  // Trapdoor dimension
  usint m_dimension;
  // gadget parameters for the trusted setup
  usint m_k0{};
  usint m_b0{};
  // gadget parameters for users
  usint m_k1{};
  usint m_b1{};
  //flag for verifying norm of signature
  bool VerifyNorm;

  /*
   *@brief Overloaded dummy method
   */
  void forceImplement() {}
};

/**
 *  @brief Class holding signing key for Ring LWE variant of GPV signing
 * algorithm with GM17 improvements. The values held in this class are trapdoor
 * and public key
 *  @tparam is the ring element
 */
template <class Element>
class GPVSignKey : public LPSignKey<Element> {
 public:
  /**
   * Default constructor
   */
  GPVSignKey() {}

  /**Constructor
   *
   * @param x trapdoor pair used for signing
   */
  explicit GPVSignKey(shared_ptr<RLWETrapdoorPair<Element>> x) {
    this->m_sk = (x);
  }

  /**
   *Destructor
   */
  ~GPVSignKey() {}

  /**
   *Method for accessing key in signing process
   *
   *@return Key used in signing
   */
  const RLWETrapdoorPair<Element>& GetSignKey() const { return *m_sk; }
  /**
   *Method for setting the private key used in the signing process
   *
   *@param &x a trapdoor pair used for signing
   */
  void SetSignKey(shared_ptr<RLWETrapdoorPair<Element>> x) { this->m_sk = (x); }

 private:
  // Trapdoor pair acting as signing key
  shared_ptr<RLWETrapdoorPair<Element>> m_sk;
  /*
   *@brief Overloaded dummy method
   */
  void forceImplement() {}
};

/**
 * @brief Class holding verification key for Ring LWE variant of GPV signing
 * algorithm with GM17 improvements. The value held in this class is the  public
 * key of the trapdoor
 * @tparam is the ring element
 */
template <class Element>
class GPVVerificationKey : public LPVerificationKey<Element> {
 public:
  /**
   *  Default constructor
   */
  GPVVerificationKey() {}

  /**
   * Constructor
   * @param vk Verification key
   */
  explicit GPVVerificationKey(shared_ptr<Matrix<Element>> vk) {
    this->m_vk = vk;
  }

  /**
   *  Destructor
   */
  ~GPVVerificationKey() {}
  /**
   *Method for accessing key in verification process
   *
   *@return Key used in verification
   */
  const Matrix<Element>& GetVerificationKey() const { return *m_vk; }
  /**
   * Method for setting key used in verification process
   *
   * @param x Key used in verification
   */
  void SetVerificationKey(shared_ptr<Matrix<Element>> x) { this->m_vk = x; }

 private:
  // Public key from trapdoor acting as verification key
  shared_ptr<Matrix<Element>> m_vk;
  /*
   *@brief Overloaded dummy method
   */
  void forceImplement() {}
};
/**
 *@brief Implementation of Ring LWE variant of GPV signature scheme. Currently
 *it supports only one type of vectors, therefore it is not templated
 *  @tparam is the ring element
 */
template <class Element>
class GPVSignatureScheme : public LPSignatureScheme<Element> {
 public:
  /**
   * Default constructor
   */
  GPVSignatureScheme() : seed(0) {}

  /**
   *Method for signing given text
   *@param m_params parameters for signing
   *@param sk private signing key
   *@param vk public verification key
   *@param pt encoding of the text to be signed
   *@param sign signature generated after the signing process - output of the
   *function
   */
  void Sign(shared_ptr<LPSignatureParameters<Element>> m_params,
            const LPSignKey<Element>& sk, const LPVerificationKey<Element>& vk,
            const LPSignPlaintext<Element>& pt, LPSignature<Element>* sign);

  void CrsGen(shared_ptr<LPSignatureParameters<Element>> m_params,
              const LPSignKey<Element>& sk, const LPVerificationKey<Element>& vk,
              const LPVerificationKey<Element>& vki, const LPSignPlaintext<Element>& usrseed,
              LPSignature<Element>* sign);
  /**
   *Method for offline perturbation sampling
   *@param m_params parameters used for signing
   *@param signKey private signing key
   *return perturbation vector
   */
  PerturbationVector<Element> SampleOffline(
      shared_ptr<LPSignatureParameters<Element>> m_params,
      const LPSignKey<Element>& signKey);

  /**
   *Method for signing given text
   *@param m_params parameters used for signing
   *@param signKey private signing key
   *@param verificationKey public verification key
   *@param Pre-computed perturbation vector
   *@param plainText encoding of the text to be signed
   *@param signatureText signature generated after the signing process - output
   *of the function
   */
  void SignOnline(shared_ptr<LPSignatureParameters<Element>> m_params,
                  const LPSignKey<Element>& signKey,
                  const LPVerificationKey<Element>& verificationKey,
                  const PerturbationVector<Element>& parturbationVector,
                  const LPSignPlaintext<Element>& pt,
                  LPSignature<Element>* signatureText);

  /**
   *Method for verifying given text & signature
   *@param m_params parameters used for the scheme
   *@param vk public verification key
   *@param sign signature to be verified
   *@param pt encoding of the text to be verified
   *@return result of the verification process
   */
  bool Verify(shared_ptr<LPSignatureParameters<Element>> m_params,
              const LPVerificationKey<Element>& vk,
              const LPSignature<Element>& sign,
              const LPSignPlaintext<Element>& pt,
              const LPSignPlaintext<Element>& seed);


  bool VerifyMulti(shared_ptr<LPSignatureParameters<Element>> m_params,
              const LPVerificationKey<Element>& vk,
              const LPSignature<Element>& sign,
              const LPSignPlaintext<Element>& pt,
              const Matrix<Element> &weight,
              const LPSignPlaintext<Element> seeds[]);

  bool VerifyAggre(shared_ptr<LPSignatureParameters<Element>> m_params,
                   const LPVerificationKey<Element> &vk,
                   const LPSignature<Element> &sign,
                   const LPSignPlaintext<Element> pts[],
                   const Matrix<Element> &weight,
                   const LPSignPlaintext<Element> seeds[]);


  void Expand_1(shared_ptr<LPSignatureParameters<Element>> sparams, const LPSignPlaintext<Element> &pt, Element &bi);
  void Expand_n(shared_ptr<LPSignatureParameters<Element>> sparams, const LPSignPlaintext<Element> &pt, Matrix<Element> &Bi);
  void Expand_hash(shared_ptr<LPSignatureParameters<Element>> sparams, const LPSignPlaintext<Element> &pt, Matrix<Element> &Bi);

  /**
   *
   *Method for generating signing and verification keys
   *@param m_params parameters used for the scheme
   *@param sk private signing key generated after trapdoor & perturbation matrix
   *- output of the function
   *@param vk public verification key generated after trapdoor - output of the
   *function
   */
  void KeyGen(shared_ptr<LPSignatureParameters<Element>> m_params,
              LPSignKey<Element>* sk, LPVerificationKey<Element>* vk, bool setup = false);

 private:
  std::vector<char> seed;
  /*
   *@brief Overloaded dummy method
   */
  void forceImplement() {}
};
}  // namespace lbcrypto
#endif
