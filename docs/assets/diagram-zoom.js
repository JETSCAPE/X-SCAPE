/* Click-to-enlarge for Mermaid diagrams.
 *
 * Material for MkDocs renders each diagram into a CLOSED shadow root attached
 * to <div class="mermaid">.  The SVG is therefore unreachable from JS — we
 * cannot query or clone it.  So instead we physically move the host element
 * (which carries its shadow root with it) into a full-screen overlay, scale it
 * up with a CSS transform, and move it back to its original spot on close.
 */
(function () {
  "use strict";

  var overlay, inner, closeBtn;
  var currentHost = null;     /* the .mermaid div currently enlarged */
  var placeholder = null;     /* marks where to put it back           */

  /* ── Build the reusable overlay once ─────────────────────────────── */
  function buildOverlay() {
    if (overlay) return;

    overlay = document.createElement("div");
    overlay.id = "dz-overlay";

    inner = document.createElement("div");
    inner.id = "dz-inner";

    closeBtn = document.createElement("button");
    closeBtn.id = "dz-close";
    closeBtn.type = "button";
    closeBtn.textContent = "✕";          /* ✕ */
    closeBtn.setAttribute("aria-label", "Close");

    overlay.appendChild(closeBtn);
    overlay.appendChild(inner);
    document.body.appendChild(overlay);

    overlay.addEventListener("click", function (e) {
      if (e.target === overlay || e.target === inner) closeZoom();
    });
    closeBtn.addEventListener("click", closeZoom);
    document.addEventListener("keydown", function (e) {
      if (e.key === "Escape") closeZoom();
    });
  }

  /* ── Enlarge a diagram ───────────────────────────────────────────── */
  function openZoom(host) {
    buildOverlay();
    if (currentHost) closeZoom();              /* safety: only one at a time */

    currentHost = host;

    /* Capture the diagram's on-page size BEFORE moving it.  In the article it
       has a *definite* width, so the inner SVG (width:100%) is laid out
       correctly.  We pin that width after moving so it can't collapse. */
    var pageRect = host.getBoundingClientRect();
    var w0 = pageRect.width;
    var h0 = pageRect.height;

    /* Remember where it lived so we can restore it exactly. */
    placeholder = document.createElement("span");
    placeholder.className = "dz-placeholder";
    placeholder.style.display = "none";
    host.parentNode.insertBefore(placeholder, host);

    /* Move the live node (shadow root included) into the overlay. */
    inner.appendChild(host);
    host.classList.add("dz-enlarged");

    /* Pin the width to its on-page value.  A width:100% SVG needs a definite
       parent width; without this it collapses to min-content (tiny) in the
       shrink-to-fit overlay.  This is the actual cause of the "smaller" bug. */
    host.style.transform = "";
    host.style.width = w0 + "px";

    /* Scale up to fill the viewport.  Based on the pinned on-page box, so a
       diagram that filled its column will fill the screen; allow <1 so a very
       tall diagram shrinks to fit rather than overflowing. */
    var scale = 1;
    if (w0 > 0 && h0 > 0) {
      scale = Math.min(
        (window.innerWidth * 0.92) / w0,
        (window.innerHeight * 0.88) / h0
      );
    }
    if (!isFinite(scale) || scale <= 0) scale = 1;
    host.style.transformOrigin = "center center";
    host.style.transform = "scale(" + scale + ")";

    overlay.classList.add("dz-visible");
    document.body.classList.add("dz-noscroll");
  }

  /* ── Restore a diagram to its original place ─────────────────────── */
  function closeZoom() {
    if (!currentHost) return;

    currentHost.classList.remove("dz-enlarged");
    currentHost.style.transform = "";
    currentHost.style.transformOrigin = "";
    currentHost.style.width = "";

    if (placeholder && placeholder.parentNode) {
      placeholder.parentNode.replaceChild(currentHost, placeholder);
    }

    placeholder = null;
    currentHost = null;

    if (overlay) overlay.classList.remove("dz-visible");
    document.body.classList.remove("dz-noscroll");
  }

  /* ── One delegated click listener handles everything ─────────────── */
  /* Clicks inside the closed shadow root are retargeted to the host in the
     light DOM, so e.target lands on (or inside) the .mermaid element.       */
  document.addEventListener("click", function (e) {
    if (overlay && overlay.contains(e.target)) return;   /* ignore overlay clicks */
    var host = e.target.closest ? e.target.closest(".mermaid") : null;
    if (host) openZoom(host);
  });

  buildOverlay();
})();
