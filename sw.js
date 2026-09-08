// Service worker: offline cache for same-origin files only. Cross-origin requests are never cached.
var CACHE = 'bottleflip-7c8388a987';
var ASSETS = ['./', './index.html', './style.css?v=7c8388a987', './guard.js?v=7c8388a987', './physics.js?v=7c8388a987', './app.js?v=7c8388a987', './manifest.webmanifest', './icon-192.png', './icon-512.png'];
self.addEventListener('install', function (e) {
  e.waitUntil(caches.open(CACHE).then(function (c) { return c.addAll(ASSETS); }).then(function () { return self.skipWaiting(); }));
});
self.addEventListener('activate', function (e) {
  e.waitUntil(caches.keys().then(function (keys) { return Promise.all(keys.filter(function (k) { return k !== CACHE; }).map(function (k) { return caches.delete(k); })); }).then(function () { return self.clients.claim(); }));
});
self.addEventListener('fetch', function (e) {
  var req = e.request;
  if (req.method !== 'GET') return;
  var url = new URL(req.url);
  if (url.origin !== self.location.origin) return;           // let the browser handle CDN/fonts with its own policies
  e.respondWith(caches.match(req, { ignoreSearch: false }).then(function (hit) {
    var fetched = fetch(req).then(function (res) {
      if (res && res.ok && res.type === 'basic') { var copy = res.clone(); caches.open(CACHE).then(function (c) { c.put(req, copy); }); }
      return res;
    }).catch(function () { return hit; });
    return hit || fetched;
  }));
});
