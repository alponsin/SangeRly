/* SangeRly - small interface helpers (no tool is changed).
   1) .sg-dlgroup: the download buttons inside become ONE menu ("Save v"). The buttons themselves are moved, not copied,
      so the downloads work exactly as before. Other elements inside (for example a format box) are kept in the menu too.
   2) a small search row (.sg-find): Enter starts the search. */
(function () {
  'use strict';

  function build(g) {
    if (g.getAttribute('data-sg-done')) return;
    var items = Array.prototype.slice.call(g.children);
    if (!items.length) return;
    g.setAttribute('data-sg-done', '1');
    var label = g.getAttribute('data-label') || 'Save';
    var ic = g.getAttribute('data-icon') || 'floppy-disk';
    g.classList.add('dropdown', 'sg-save');
    g.removeAttribute('style');

    var btn = document.createElement('button');
    btn.type = 'button';
    btn.className = 'btn btn-sm btn-outline-secondary dropdown-toggle sg-save-btn';
    btn.setAttribute('data-bs-toggle', 'dropdown');
    btn.setAttribute('data-bs-auto-close', 'outside');
    btn.setAttribute('aria-expanded', 'false');
    btn.innerHTML = '<i class="fa fa-' + ic + '" role="presentation"></i> <span class="sg-save-lbl"></span>';
    btn.querySelector('.sg-save-lbl').textContent = label;

    var ul = document.createElement('ul');
    ul.className = 'dropdown-menu sg-save-menu';
    items.forEach(function (el) {
      var li = document.createElement('li');
      var isLink = el.matches && el.matches('a, button');
      if (isLink) {
        el.className = 'dropdown-item ' + (el.classList.contains('shiny-download-link') ? 'shiny-download-link' : '') +
                       (el.classList.contains('action-button') ? ' action-button' : '');
        li.appendChild(el);
      } else {
        li.className = 'dropdown-item-text';
        li.appendChild(el);
      }
      ul.appendChild(li);
    });
    g.innerHTML = '';
    g.appendChild(btn);
    g.appendChild(ul);

    try {   // the menu may sit inside a box that scrolls: let it float above everything
      if (window.bootstrap && window.bootstrap.Dropdown) {
        new window.bootstrap.Dropdown(btn, { popperConfig: function (d) { d.strategy = 'fixed'; return d; } });
      }
    } catch (e) { /* the plain menu still works */ }

    ul.addEventListener('click', function (e) {   // after a download or an action the menu closes
      var a = e.target.closest ? e.target.closest('a.dropdown-item, button.dropdown-item') : null;
      if (!a) return;
      setTimeout(function () {
        try { window.bootstrap.Dropdown.getOrCreateInstance(btn).hide(); } catch (err) { /* ignore */ }
      }, 150);
    });
  }

  var pending = false;
  function scan() {
    pending = false;
    var gs = document.querySelectorAll('.sg-dlgroup:not([data-sg-done])');
    for (var i = 0; i < gs.length; i++) build(gs[i]);
  }
  function schedule() {
    if (pending) return;
    pending = true;
    (window.requestAnimationFrame || window.setTimeout)(scan, 0);
  }

  function start() {
    scan();
    new MutationObserver(schedule).observe(document.documentElement, { childList: true, subtree: true });
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();

  // Enter in the search row = press the main button of that row
  document.addEventListener('keydown', function (e) {
    if (e.key !== 'Enter' || !e.target || !e.target.matches || !e.target.matches('.sg-find input[type=text]')) return;
    var row = e.target.closest('.sg-find');
    var b = row && row.querySelector('button.btn-primary');
    if (!b) return;
    e.preventDefault();
    if (window.jQuery) window.jQuery(e.target).trigger('change');   // send the typed text first
    setTimeout(function () { b.click(); }, 80);
  });
})();
