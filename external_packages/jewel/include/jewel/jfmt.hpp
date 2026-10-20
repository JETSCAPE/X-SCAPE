// JEWEL 2.6.0 C++ port -- std::format / std::println, or a small replacement.
//
// The port uses std::format and std::println (C++20/23). Inside JETSCAPE it is
// also built with older standard libraries (e.g. GCC 11, which has neither), so
// all formatting goes through jewel::jfmt. With a library that provides
// <format> this is std::format itself. Otherwise a minimal formatter is used
// that supports what the port needs: automatic "{}" fields with an optional
// spec [[fill]align][0][width][.precision][d|e|f|g|s] for integers, floating
// point values, characters, booleans and strings, plus "{{" / "}}" escapes.
// Default floating-point output is the shortest round-trip form, as in
// std::format.
#pragma once

#include <cstdio>
#include <string>
#include <string_view>
#include <version>

#if defined(__cpp_lib_format) && !defined(JEWEL_JFMT_FORCE_SHIM)
#include <format>
#if defined(__cpp_lib_print)
#include <print>
#endif

namespace jewel::jfmt {
using std::format;
using std::format_string;
#if defined(__cpp_lib_print)
using std::println;
#else
template <typename... Args>
void println(std::format_string<Args...> f, Args&&... args) {
  std::string s = std::format(f, std::forward<Args>(args)...);
  s += '\n';
  std::fputs(s.c_str(), stdout);
}
template <typename... Args>
void println(std::FILE* out, std::format_string<Args...> f, Args&&... args) {
  std::string s = std::format(f, std::forward<Args>(args)...);
  s += '\n';
  std::fputs(s.c_str(), out);
}
#endif
}  // namespace jewel::jfmt

#else  // ---- replacement formatter -------------------------------------------

#include <charconv>
#include <stdexcept>
#include <type_traits>

namespace jewel::jfmt {

// Format strings are not checked at compile time here.
template <typename... Args>
using format_string = std::string_view;

namespace detail {

struct Spec {
  char fill = ' ';
  char align = 0;  // '<', '>', '^' or 0 (type default)
  bool zero = false;
  int width = -1;
  int prec = -1;
  char type = 0;
};

inline Spec parse_spec(std::string_view s) {
  Spec sp;
  std::size_t i = 0;
  auto is_align = [](char c) { return c == '<' || c == '>' || c == '^'; };
  if (s.size() >= 2 && is_align(s[1])) {
    sp.fill = s[0];
    sp.align = s[1];
    i = 2;
  } else if (!s.empty() && is_align(s[0])) {
    sp.align = s[0];
    i = 1;
  }
  if (i < s.size() && s[i] == '0') {
    sp.zero = true;
    ++i;
  }
  int w = 0;
  bool hasw = false;
  while (i < s.size() && s[i] >= '0' && s[i] <= '9') {
    w = 10 * w + (s[i++] - '0');
    hasw = true;
  }
  if (hasw) sp.width = w;
  if (i < s.size() && s[i] == '.') {
    ++i;
    int p = 0;
    while (i < s.size() && s[i] >= '0' && s[i] <= '9') p = 10 * p + (s[i++] - '0');
    sp.prec = p;
  }
  if (i < s.size()) sp.type = s[i++];
  if (i != s.size()) throw std::runtime_error("jfmt: unsupported format spec");
  return sp;
}

inline void pad(std::string& out, std::string body, const Spec& sp, bool numeric) {
  const int n = static_cast<int>(body.size());
  if (sp.width <= n) {
    out += body;
    return;
  }
  const int fill = sp.width - n;
  if (numeric && sp.zero && sp.align == 0) {
    // zero padding goes after the sign
    std::size_t sign = (!body.empty() && (body[0] == '-' || body[0] == '+')) ? 1 : 0;
    body.insert(sign, static_cast<std::size_t>(fill), '0');
    out += body;
    return;
  }
  const char align = sp.align ? sp.align : (numeric ? '>' : '<');
  if (align == '<') {
    out += body;
    out.append(static_cast<std::size_t>(fill), sp.fill);
  } else if (align == '>') {
    out.append(static_cast<std::size_t>(fill), sp.fill);
    out += body;
  } else {
    out.append(static_cast<std::size_t>(fill / 2), sp.fill);
    out += body;
    out.append(static_cast<std::size_t>(fill - fill / 2), sp.fill);
  }
}

template <typename T>
std::string float_body(T v, const Spec& sp) {
  char buf[512];
  if (sp.prec < 0 && sp.type == 0) {
    auto r = std::to_chars(buf, buf + sizeof(buf), v);
    return std::string(buf, r.ptr);
  }
  const int prec = sp.prec < 0 ? 6 : sp.prec;
  const char* f = "%.*g";
  if (sp.type == 'f' || sp.type == 'F') f = "%.*f";
  else if (sp.type == 'e') f = "%.*e";
  else if (sp.type == 'E') f = "%.*E";
  std::snprintf(buf, sizeof(buf), f, prec, static_cast<double>(v));
  return buf;
}

template <typename T>
void format_one(std::string& out, const T& v, const Spec& sp) {
  using U = std::remove_cv_t<std::remove_reference_t<T>>;
  if constexpr (std::is_same_v<U, bool>) {
    pad(out, v ? "true" : "false", sp, false);
  } else if constexpr (std::is_same_v<U, char>) {
    pad(out, std::string(1, v), sp, false);
  } else if constexpr (std::is_integral_v<U>) {
    pad(out, std::to_string(v), sp, true);
  } else if constexpr (std::is_floating_point_v<U>) {
    pad(out, float_body(v, sp), sp, true);
  } else {
    pad(out, std::string(std::string_view(v)), sp, false);
  }
}

template <typename T>
void format_nth(std::string& out, std::size_t n, const Spec& sp, const T& first) {
  if (n != 0) throw std::runtime_error("jfmt: too few arguments");
  format_one(out, first, sp);
}
template <typename T, typename... Rest>
void format_nth(std::string& out, std::size_t n, const Spec& sp, const T& first,
                const Rest&... rest) {
  if (n == 0) format_one(out, first, sp);
  else format_nth(out, n - 1, sp, rest...);
}

}  // namespace detail

template <typename... Args>
std::string format(std::string_view f, const Args&... args) {
  std::string out;
  out.reserve(f.size() + 16 * sizeof...(Args));
  std::size_t next = 0;
  for (std::size_t i = 0; i < f.size(); ++i) {
    const char c = f[i];
    if (c == '{') {
      if (i + 1 < f.size() && f[i + 1] == '{') {
        out += '{';
        ++i;
        continue;
      }
      const std::size_t close = f.find('}', i);
      if (close == std::string_view::npos) throw std::runtime_error("jfmt: unmatched '{'");
      std::string_view field = f.substr(i + 1, close - i - 1);
      const std::size_t colon = field.find(':');
      if (colon != 0 && !field.empty() && colon != std::string_view::npos)
        throw std::runtime_error("jfmt: explicit argument indices are not supported");
      const detail::Spec sp =
          colon == std::string_view::npos ? detail::Spec{} : detail::parse_spec(field.substr(colon + 1));
      if constexpr (sizeof...(Args) > 0) {
        detail::format_nth(out, next++, sp, args...);
      } else {
        throw std::runtime_error("jfmt: too few arguments");
      }
      i = close;
    } else if (c == '}') {
      if (i + 1 < f.size() && f[i + 1] == '}') ++i;
      out += '}';
    } else {
      out += c;
    }
  }
  return out;
}

template <typename... Args>
void println(std::string_view f, const Args&... args) {
  std::string s = format(f, args...);
  s += '\n';
  std::fputs(s.c_str(), stdout);
}
template <typename... Args>
void println(std::FILE* out, std::string_view f, const Args&... args) {
  std::string s = format(f, args...);
  s += '\n';
  std::fputs(s.c_str(), out);
}

}  // namespace jewel::jfmt

#endif
