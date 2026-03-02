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

// Jetscape final state {hadrons,partons} writer ascii class
// Based on JetScapeWriterStream.
// author: Raymond Ehlers <raymond.ehlers@cern.ch>, LBL/UCB

#ifndef JETSCAPEWRITERFINALSTATESTREAM_H
#define JETSCAPEWRITERFINALSTATESTREAM_H

#include <fstream>
#include <string>

#ifdef USE_GZIP
#include "gzstream.h"
#endif

#include "JetScapeWriter.h"

using std::ofstream;

namespace Jetscape {

/**
 * @class JetScapeWriterFinalStateStream
 * @brief Template base class for writing final state particles to streams.
 *
 * This class writes either hadrons or partons to an ASCII (or gzipped) stream.
 * The specialization (hadrons or partons) is determined in derived classes.
 *
 * @tparam T Output stream type (e.g. std::ofstream or ogzstream).
 */
template <class T>
class JetScapeWriterFinalStateStream : public JetScapeWriter {
 public:
  /**
   * @brief Default constructor.
   */
  JetScapeWriterFinalStateStream<T>(){};

  /**
   * @brief Constructor with output file name.
   * @param m_file_name_out Name of the output file.
   */
  JetScapeWriterFinalStateStream<T>(string m_file_name_out);

  /**
   * @brief Destructor.
   */
  virtual ~JetScapeWriterFinalStateStream<T>();

  /**
   * @brief Initialize the writer.
   */
  void InitTask();

  /**
   * @brief Execute writing routine.
   */
  void ExecuteTask();

  /**
   * @brief Get the name of the writer type.
   * @throws std::runtime_error if called on the base class.
   */
  virtual std::string GetName() {
    throw std::runtime_error("Don't use the base class");
  }

  /**
   * @brief Check the status of the output file.
   * @return True if the file is ready for writing.
   */
  bool GetStatus() { return output_file.good(); }

  /**
   * @brief Close the output file.
   *
   * This function is also used to add cross-section and error information.
   */
  void Close();

  /**
   * @brief Write a parton shower to the stream.
   * @param ps Weak pointer to a PartonShower object.
   */
  void Write(weak_ptr<PartonShower> ps);

  /**
   * @brief Write a hadron to the stream.
   * @param h Weak pointer to a Hadron object.
   */
  void Write(weak_ptr<Hadron> h);
  // We aren't interested in the individual partons or vertices, so skip them.

  /**
   * @brief Write header information to the stream.
   *
   * No-op in this implementation.
   */
  void WriteHeaderToFile(){};

  /**
   * @brief Write event-level information to the stream.
   */
  void WriteEvent();

  /**
   * @brief Write a raw string to the stream.
   * @param s String to be written.
   */
  void Write(string s) { output_file << s << endl; }
  // Intentionally make these no-ops since we want to fully control our output
  // from this class. Tasks will often directly call these functions, so we need
  // to prevent them from doing so.

  /**
   * @brief Write a comment to the stream.
   *
   * No-op in this implementation to enforce controlled formatting.
   * @param s Comment string (ignored).
   */
  void WriteComment(string s) {}

  /**
   * @brief Write whitespace to the stream.
   *
   * No-op in this implementation to enforce controlled formatting.
   * @param s Whitespace string (ignored).
   */
  void WriteWhiteSpace(string s) {}

 protected:
  T output_file;               ///< Output file stream.
  unsigned int headerVersion;  ///< Version number for the header format.

  /** Collected particles to write. */
  std::vector<std::shared_ptr<JetScapeParticleBase>> particles;
  bool writeCentrality;  ///< Whether to write centrality information.
  bool writePtHat;       ///< Whether to write pTHat information.

  /** Particle status codes to skip when writing. */
  std::vector<int> particleStatusToSkip;
  bool binaryOutput;  ///< Whether to write in binary format.
};

/**
 * @class JetScapeWriterFinalStatePartonsStream
 * @brief Specialization for writing final state partons to a stream.
 *
 * This class inherits from JetScapeWriterFinalStateStream and provides
 * functionality specific to parton output. Hadrons are explicitly skipped.
 *
 * @tparam T Output stream type.
 */
template <class T>
class JetScapeWriterFinalStatePartonsStream
    : public JetScapeWriterFinalStateStream<T> {
  /**
   * @brief Get the name of this writer type.
   * @return The string "partons".
   */
  std::string GetName() { return "partons"; }

  /**
   * @brief Skip writing hadrons.
   * @param h Weak pointer to a Hadron (ignored).
   * @note Don't collect the hadrons by making this a no-op
   */
  void Write(weak_ptr<Hadron> h) {}

 protected:
  /**
   * @brief Register the parton writer with the JETSCAPE framework.
   */
  static RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ofstream>>
      regParton;

  /**
   * @brief Register the gzipped parton writer with the JETSCAPE framework.
   */
  static RegisterJetScapeModule<
      JetScapeWriterFinalStatePartonsStream<ogzstream>>
      regPartonGZ;
};

/**
 * @class JetScapeWriterFinalStateHadronsStream
 * @brief Specialization for writing final state hadrons to a stream.
 *
 * This class inherits from JetScapeWriterFinalStateStream and provides
 * functionality specific to hadron output. Parton showers are explicitly
 * skipped.
 *
 * @tparam T Output stream type.
 */
template <class T>
class JetScapeWriterFinalStateHadronsStream
    : public JetScapeWriterFinalStateStream<T> {
  /**
   * @brief Get the name of this writer type.
   * @return The string "hadrons".
   */
  std::string GetName() { return "hadrons"; }

  /**
   * @brief Skip writing parton showers.
   * @param ps Weak pointer to a PartonShower (ignored).
   * @note Don't collect the parton showers by making this a no-op
   */
  void Write(weak_ptr<PartonShower> ps) {}

 protected:
  /**
   * @brief Register the hadron writer with the JETSCAPE framework.
   */
  static RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ofstream>>
      regHadron;
  /**
   * @brief Register the gzipped hadron writer with the JETSCAPE framework.
   */
  static RegisterJetScapeModule<
      JetScapeWriterFinalStateHadronsStream<ogzstream>>
      regHadronGZ;
};

/**
 * @typedef JetScapeWriterFinalStateStreamAscii
 * @brief Alias for the ASCII final state writer.
 */
typedef JetScapeWriterFinalStateStream<ofstream>
    JetScapeWriterFinalStateStreamAscii;

/**
 * @typedef JetScapeWriterFinalStatePartonsAscii
 * @brief Alias for the ASCII parton writer.
 */
typedef JetScapeWriterFinalStatePartonsStream<ofstream>
    JetScapeWriterFinalStatePartonsAscii;

/** @typedef JetScapeWriterFinalStateHadronsAscii
 * @brief Alias for the ASCII hadron writer.
 */
typedef JetScapeWriterFinalStateHadronsStream<ofstream>
    JetScapeWriterFinalStateHadronsAscii;
#ifdef USE_GZIP

/** @typedef JetScapeWriterFinalStateStreamAsciiGZ
 * @brief Alias for the gzipped ASCII final state writer.
 */
typedef JetScapeWriterFinalStateStream<ogzstream>
    JetScapeWriterFinalStateStreamAsciiGZ;

/** @typedef JetScapeWriterFinalStatePartonsAsciiGZ
 * @brief Alias for the gzipped ASCII parton writer.
 */
typedef JetScapeWriterFinalStatePartonsStream<ogzstream>
    JetScapeWriterFinalStatePartonsAsciiGZ;

/** @typedef JetScapeWriterFinalStateHadronsAsciiGZ
 * @brief Alias for the gzipped ASCII hadron writer.
 */
typedef JetScapeWriterFinalStateHadronsStream<ogzstream>
    JetScapeWriterFinalStateHadronsAsciiGZ;
#endif

}  // end namespace Jetscape

#endif  // JETSCAPEWRITERFINALSTATESTREAM_H
