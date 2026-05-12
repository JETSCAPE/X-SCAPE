/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

// Jetscape Qn vector writer ascii class
// Based on JetScapeWriterStream, JetScapeWriterFinalStateStream
// author: Xiang-Yu Wu <xiangyu.wu2@mail.mcgill.ca>

#ifndef JETSCAPEWRITERSTREAM_H
#define JETSCAPEWRITERSTREAM_H

#include <fstream>
#include <string>

#ifdef USE_GZIP
#include "gzstream.h"
#endif

#include "JetScapeWriter.h"

using std::ofstream;

#include "pdgcode.h"

namespace Jetscape {

/**
 * @class JetScapeWriterQnVectorStream
 * @brief Writer module that accumulates hadrons and computes Qn-vector
 * observables.
 *
 * This templated writer collects all final-state hadrons from an event,
 * bins them in transverse momentum and rapidity, and computes flow
 * coefficients (Qn-vectors) for specified particle species or charge classes.
 *
 * The results are written in a structured tabular ASCII (or compressed)
 * format that includes per-event metadata and statistical errors.
 *
 * The configuration (pT range, number of bins, rapidity range, harmonic
 * order) is read from the XML input file (e.g. `QnVector_pTmin`,
 * `QnVector_NpT`, `QnVector_Norder`, etc.).
 *
 * @tparam T Output stream type, typically `std::ofstream` or `ogzstream`.
 */

template <class T>
class JetScapeWriterQnVectorStream : public JetScapeWriter {
 public:
  /**
   * @brief Default constructor.
   */
  JetScapeWriterQnVectorStream<T>(){};

  /**
   * @brief Construct a writer with an explicit output file name.
   * @param m_file_name_out Name of the file to be opened for writing.
   */
  JetScapeWriterQnVectorStream<T>(string m_file_name_out);

  /**
   * @brief Destructor.
   *
   * Ensures the output file is closed and writes an event footer if the
   * module is still active.
   */
  virtual ~JetScapeWriterQnVectorStream<T>();

  /**
   * @brief Initialize the writer.
   *
   * Reads binning and Qn-vector configuration from the XML input,
   * opens the output file, and writes a file-level header.
   */
  void Init();

  /**
   * @brief Execute the module.
   *
   * Currently a placeholder; event writing is controlled by
   * `WriteEvent()`.
   */
  void Exec();

  /**
   * @brief Return the name of the module.
   * @return Always returns `"QnVector"`.
   */
  virtual std::string GetName() { return "QnVector"; }

  /**
   * @brief Query whether the output file is in a good state.
   * @return `true` if the file stream is valid, `false` otherwise.
   */
  bool GetStatus() { return output_file.good(); }

  /**
   * @brief Close the output stream.
   *
   * Writes an event footer marking the end of the file and closes
   * the underlying stream.
   */
  void Close();

  /**
   * @brief Write a parton shower (no-op for this writer).
   */
  void Write(weak_ptr<PartonShower> ps){};

  /**
   * @brief Add a hadron to the buffer for Qn-vector analysis.
   * @param h Weak pointer to a hadron.
   */
  void Write(weak_ptr<Hadron> h);
  // We aren't interested in the individual partons or vertices, so skip them.

  /**
   * @brief Write a file header (no-op, handled by Init()).
   */
  void WriteHeaderToFile(){};

  /**
   * @brief Process all accumulated hadrons and write Qn-vectors to file.
   *
   * Loops over selected particle species, fills Qn-vectors using
   * `Qvector` objects and writes results in tabular form.
   */
  void WriteEvent();

  /**
   * @brief Write a raw string to the output stream.
   * @param s String to write.
   */
  void Write(string s) { output_file << s << endl; }

  // Intentionally make these no-ops since we want to fully control our output
  // from this class. Tasks will often directly call these functions, so we need
  // to prevent them from doing so.

  /**
   * @brief Suppress comment writing (no-op).
   *
   * Overridden to prevent tasks from writing arbitrary comments.
   */
  void WriteComment(string s) {}

  /**
   * @brief Suppress whitespace writing (no-op).
   */
  void WriteWhiteSpace(string s) {}

 protected:
  T output_file;  //!< Output file stream (ASCII or GZip).
  std::vector<std::shared_ptr<Hadron>>
      particles;         //!< Collected hadrons for the current event.
  bool writeCentrality;  //!< Whether to write event centrality.
  static RegisterJetScapeModule<JetScapeWriterQnVectorStream<ofstream>>
      regQnVector;  /// Registration of the GZip writer module (if enabled).
  static RegisterJetScapeModule<JetScapeWriterQnVectorStream<ogzstream>>
      regQnVectorGZ;

 private:
  double pTmin_;              //!< Minimum pT for binning.
  double pTmax_;              //!< Maximum pT for binning.
  double rapmin_;             //!< Minimum rapidity for binning.
  double rapmax_;             //!< Maximum rapidity for binning.
  int npT_;                   //!< Number of pT bins.
  int nrap_;                  //!< Number of rapidity bins.
  int norder_;                //!< Maximum harmonic order for Qn-vectors.
  std::map<int, int> chpdg_;  //!< PDG code lookup for charged hadrons.
};

/**
 * @typedef JetScapeWriterQnVectorAscii
 * @brief Convenience typedef for ASCII Qn-vector writer (`ofstream` backend).
 */
typedef JetScapeWriterQnVectorStream<ofstream> JetScapeWriterQnVectorAscii;

#ifdef USE_GZIP
/**
 * @typedef JetScapeWriterQnVectorAsciiGZ
 * @brief Convenience typedef for GZip-compressed Qn-vector writer (`ogzstream`
 * backend).
 */
typedef JetScapeWriterQnVectorStream<ogzstream> JetScapeWriterQnVectorAsciiGZ;
#endif

}  // end namespace Jetscape

#endif  // JETSCAPEWRITERSTREAM_H
