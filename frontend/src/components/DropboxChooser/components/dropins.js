/* eslint-disable */
export default function dropins() {
  (function (e, o) {
    for (const t in o) e[t] = o[t];
  })(
    window,
    (function (e) {
      const o = {};
      function t(n) {
        if (o[n]) return o[n].exports;
        const r = (o[n] = { i: n, l: !1, exports: {} });
        return e[n].call(r.exports, r, r.exports, t), (r.l = !0), r.exports;
      }
      return (
        (t.m = e),
        (t.c = o),
        (t.d = function (e, o, n) {
          t.o(e, o) || Object.defineProperty(e, o, { enumerable: !0, get: n });
        }),
        (t.r = function (e) {
          "undefined" != typeof Symbol &&
            Symbol.toStringTag &&
            Object.defineProperty(e, Symbol.toStringTag, { value: "Module" }),
            Object.defineProperty(e, "__esModule", { value: !0 });
        }),
        (t.t = function (e, o) {
          if ((1 & o && (e = t(e)), 8 & o)) return e;
          if (4 & o && "object" == typeof e && e && e.__esModule) return e;
          const n = Object.create(null);
          if (
            (t.r(n),
            Object.defineProperty(n, "default", { enumerable: !0, value: e }),
            2 & o && "string" != typeof e)
          )
            for (const r in e)
              t.d(
                n,
                r,
                function (o) {
                  return e[o];
                }.bind(null, r)
              );
          return n;
        }),
        (t.n = function (e) {
          const o =
            e && e.__esModule
              ? function () {
                  return e.default;
                }
              : function () {
                  return e;
                };
          return t.d(o, "a", o), o;
        }),
        (t.o = function (e, o) {
          return Object.prototype.hasOwnProperty.call(e, o);
        }),
        (t.p = ""),
        t((t.s = 6))
      );
    })([
      function (e, o, t) {
        var n,
          r =
            (this && this.__assign) ||
            function () {
              return (r =
                Object.assign ||
                function (e) {
                  for (var o, t = 1, n = arguments.length; t < n; t++)
                    for (const r in (o = arguments[t]))
                      Object.prototype.hasOwnProperty.call(o, r) &&
                        (e[r] = o[r]);
                  return e;
                }).apply(this, arguments);
            };
        function i(e, o, t) {
          return "" + e + (-1 === e.indexOf("?") ? "?" : "&") + o + "=" + t;
        }
        function s(e) {
          return i(e, "version", encodeURIComponent(Dropbox.VERSION));
        }
        function a(e, o) {
          const t = encodeURIComponent(
              window.location.protocol + "//" + window.location.host
            ),
            n = encodeURIComponent(Dropbox.appKey),
            r = encodeURIComponent(e.linkType || ""),
            a = encodeURIComponent(e._trigger || "js"),
            l = Boolean(e.multiselect),
            c = encodeURIComponent(
              d(e.extensions, "join", function (e) {
                return e.join(" ");
              }) || ""
            ),
            u = Boolean(e.folderselect);
          o = Boolean(o);
          let p =
            Dropbox.baseUrl +
            "/chooser?origin=" +
            t +
            "&app_key=" +
            n +
            "&link_type=" +
            r +
            "&trigger=" +
            a +
            "&multiselect=" +
            l +
            "&extensions=" +
            c +
            "&folderselect=" +
            u +
            "&iframe=" +
            o;
          return (
            void 0 !== e.fileselect &&
              (p = i(p, "fileselect", Boolean(e.fileselect))),
            void 0 !== e.sizeLimit && (p = i(p, "size_limit", e.sizeLimit)),
            null != e.initialNavigation &&
              (null != e.initialNavigation.mode &&
                (p = i(
                  p,
                  "initial_navigation_mode",
                  encodeURIComponent(e.initialNavigation.mode)
                )),
              null != e.initialNavigation.role &&
                (p = i(
                  p,
                  "initial_navigation_role",
                  encodeURIComponent(e.initialNavigation.role)
                )),
              e.initialNavigation.cursor &&
                (p = i(
                  p,
                  "initial_navigation_cursor",
                  encodeURIComponent(e.initialNavigation.cursor)
                ))),
            null != e.initialViewType &&
              (p = i(
                p,
                "initial_view_type",
                encodeURIComponent(e.initialViewType)
              )),
            null != e.fields &&
              (p = i(
                p,
                "fields",
                encodeURIComponent(
                  "function" == typeof e.fields.join
                    ? e.fields.join(" ")
                    : void 0
                )
              )),
            !1 === e.showSignOut && (p = i(p, "show_sign_out", "false")),
            e.initialNavigationPath &&
              (p = i(
                p,
                "initial_navigation_path",
                encodeURIComponent(e.initialNavigationPath)
              )),
            e.requiredPermissions &&
              (p = i(
                p,
                "required_permissions",
                encodeURIComponent(
                  d(e.requiredPermissions, "join", function (e) {
                    return e.join(" ");
                  }) || ""
                )
              )),
            s(p)
          );
        }
        function l(e) {
          var t = {
            options: r(r({}, e), {
              success: function (n, r) {
                "function" == typeof e.success && e.success(n, r),
                  o.currentChooserSession === t &&
                    (o.currentChooserSession = null);
              },
              cancel: function (n) {
                "function" == typeof e.cancel && e.cancel(n),
                  o.currentChooserSession === t &&
                    (o.currentChooserSession = null);
              },
            }),
          };
          return (o.currentChooserSession = t), t;
        }
        function c(e) {
          const o = document.createElement("iframe");
          return (
            (o.src = "about:blank"),
            (o._postAction = e),
            (o.name = "dropbox-dropins"),
            (o.style.display = "block"),
            (o.style.backgroundColor = "white"),
            (o.style.border = "none"),
            o
          );
        }
        function u(e, o) {
          let t,
            r = encodeURIComponent(Dropbox.appKey),
            i =
              Dropbox.baseUrl +
              "/dropins/job_status?job=" +
              o +
              "&app_key=" +
              r;
          i = s(i);
          const a = function (o) {
            if ("COMPLETE" === o.status) {
              if (
                ("function" == typeof e.progress && e.progress(1),
                "function" == typeof e.success)
              )
                if (0 !== e.success.length && o.file_ids) {
                  const n = { fileIds: o.file_ids };
                  e.success(n);
                } else e.success();
            } else
              "PENDING" === o.status || "DOWNLOADING" === o.status
                ? (null != o.progress &&
                    "function" == typeof e.progress &&
                    e.progress(o.progress / 100),
                  setTimeout(t, 1500))
                : "FAILED" === o.status &&
                  "function" == typeof e.error &&
                  e.error(o.error);
          };
          if ("withCredentials" in new XMLHttpRequest())
            t = function () {
              const o = new XMLHttpRequest();
              return (
                (o.onload = function () {
                  return a(JSON.parse(o.responseText));
                }),
                (o.onerror = function () {
                  return "function" == typeof e.error ? e.error() : void 0;
                }),
                o.open("GET", i, !0),
                o.send()
              );
            };
          else if (Dropbox.disableJSONP) {
            if (
              "undefined" == typeof XDomainRequest ||
              null === XDomainRequest ||
              "https:" !== document.location.protocol
            )
              throw new Error(
                "Unable to find suitable means of cross domain communication"
              );
            t = function () {
              const o = new XDomainRequest();
              return (
                (o.onload = function () {
                  return a(JSON.parse(o.responseText));
                }),
                (o.onerror = function () {
                  return "function" == typeof e.error ? e.error() : void 0;
                }),
                o.open("get", i),
                o.send()
              );
            };
          } else
            t = function () {
              let o = "DropboxJsonpCallback" + n++,
                t = !1;
              window[o] = function (e) {
                return (t = !0), a(e);
              };
              const r = document.createElement("script");
              return (
                (r.src = i + "&callback=" + o),
                (r.onreadystatechange = function () {
                  if ("loaded" === r.readyState)
                    return (
                      t || ("function" == typeof e.error && e.error()),
                      null != r.parentNode
                        ? r.parentNode.removeChild(r)
                        : void 0
                    );
                }),
                document.getElementsByTagName("head")[0].appendChild(r)
              );
            };
          return "function" == typeof e.progress && e.progress(0), t();
        }
        function p(e, t, n) {
          let r,
            i = JSON.parse(e.data);
          switch (
            ((r =
              null != o.ieframe && n._popup
                ? o.ieframe.contentWindow
                : e.source),
            void 0 !== i.sequence_number &&
              r.postMessage(
                JSON.stringify({
                  method: "ack",
                  sequence_number: i.sequence_number,
                }),
                Dropbox.baseUrl
              ),
            i.method)
          ) {
            case "origin_request":
              e.source.postMessage(
                JSON.stringify({ method: "origin" }),
                Dropbox.baseUrl
              );
              break;
            case "ready":
              if (null != n.files) {
                let s = void 0;
                if (n._fetch_url_on_save) {
                  for (var a = [], l = 0; l < n.files.length; l++) {
                    const c = n.files[l];
                    a.push({ filename: c.filename });
                  }
                  s = JSON.stringify({
                    method: "files_with_callback",
                    params: a,
                  });
                } else s = JSON.stringify({ method: "files", params: n.files });
                if (
                  (r.postMessage(s, Dropbox.baseUrl), null != n._ews_auth_token)
                ) {
                  const p = JSON.stringify({
                    method: "ews_auth_token",
                    params: { ews_auth_token: n._ews_auth_token },
                  });
                  r.postMessage(p, Dropbox.baseUrl);
                }
              }
              "function" == typeof n.ready && n.ready();
              break;
            case "files_selected":
            case "files_saved":
              "function" == typeof t && t(),
                "function" == typeof n.success &&
                  n.success(i.params, o.last_navigation);
              break;
            case "cursor_changed":
              o.last_navigation = { cursor: i.params };
              break;
            case "progress":
              "function" == typeof n.progress && n.progress(i.params);
              break;
            case "close_dialog":
              "function" == typeof t && t(),
                "function" == typeof n.cancel && n.cancel(o.last_navigation);
              break;
            case "resize":
              "function" == typeof n.resize && n.resize(i.params);
              break;
            case "error":
              "function" == typeof n.error && n.error(i.params);
              break;
            case "error_and_close":
              "function" == typeof t && t(),
                "function" == typeof n.error && n.error(i.params);
              break;
            case "job_id":
              "function" == typeof t && t(), u(n, i.params);
              break;
            case "save_callback":
              !(function (e, o, t) {
                if (e._fetch_url_on_save) {
                  const n = e.fetch_urls_fn;
                  "function" != typeof n &&
                    "function" == typeof e.error &&
                    e.error(
                      "Something went wrong, file url callback not provided."
                    ),
                    n(function (e) {
                      if (null == e)
                        throw new Error(
                          "Please supply {urls: [...]} to success callback"
                        );
                      if (null != e.url && null != e.urls)
                        throw new Error(
                          "Do not pass both url and urls to the save callback"
                        );
                      if ((null != e.url && (e.urls = [e.url]), null == e.urls))
                        throw new Error(
                          "Please supply {urls: [...]} to success callback"
                        );
                      return (
                        (i = {
                          method: "continue_saving",
                          params: { download_urls: e.urls },
                        }),
                        void r.postMessage(JSON.stringify(i), Dropbox.baseUrl)
                      );
                    }, o);
                }
              })(n, i.params);
              break;
            case "_debug_log":
              "undefined" != typeof console &&
                null !== console &&
                console.log(i.params.msg);
          }
        }
        function d(e, o, t) {
          return null != e && "function" == typeof e[o] ? t(e, o) : void 0;
        }
        Object.defineProperty(o, "__esModule", { value: !0 }),
          null == window.Dropbox && (window.Dropbox = {}),
          (o.popupDimensionsString = function (e, o) {
            return (
              "width=" +
              e +
              ",height=" +
              o +
              ",left=" +
              (window.screenX +
                ((window.outerWidth || document.documentElement.offsetWidth) -
                  e) /
                  2) +
              ",top=" +
              (window.screenY +
                ((window.outerHeight || document.documentElement.offsetHeight) -
                  o) /
                  2)
            );
          }),
          (o.chooserUrl = a),
          (o.createIEFrame = function () {
            /\bTrident\b/.test(navigator.userAgent) &&
              null != document.body &&
              null == o.ieframe &&
              ((o.ieframe = document.createElement("iframe")),
              o.ieframe.setAttribute("id", "dropbox_xcomm"),
              o.ieframe.setAttribute(
                "src",
                Dropbox.baseUrl + "/static/api/1/xcomm.html"
              ),
              (o.ieframe.style.display = "none"),
              document.body.appendChild(o.ieframe));
          }),
          (o.createChooserSession = l),
          (o.createWidgetElement = c),
          (o.handleJobId = u),
          (o.handleMessageEvent = p),
          (o.saverUrl = function (e) {
            const o = encodeURIComponent(
                window.location.protocol + "//" + window.location.host
              ),
              t = encodeURIComponent(Dropbox.appKey);
            return (
              (e = Boolean(e)),
              s(
                Dropbox.baseUrl +
                  "/saver?origin=" +
                  o +
                  "&app_key=" +
                  t +
                  "&iframe=" +
                  e
              )
            );
          }),
          (o.__guardMethod__ = d),
          (o.initModule = function () {
            (o.last_navigation = {}),
              (o.ieframe = null),
              (o.currentChooserSession = null),
              (n = 1),
              null == Dropbox.baseUrl &&
                (Dropbox.baseUrl = "https://www.dropbox.com"),
              null == Dropbox.blockBaseUrl &&
                (Dropbox.blockBaseUrl = "https://dl-web.dropbox.com"),
              (Dropbox.addListener = function (e, o, t) {
                e.addEventListener
                  ? e.addEventListener(o, t, !1)
                  : e.attachEvent("on" + o, function (e) {
                      return (
                        (e.preventDefault = function () {
                          return !1;
                        }),
                        t(e)
                      );
                    });
              }),
              (Dropbox.removeListener = function (e, o, t) {
                e.removeEventListener
                  ? e.removeEventListener(o, t, !1)
                  : e.detachEvent("on" + o, t);
              }),
              (Dropbox.createChooserWidget = function (e) {
                const o = l(e),
                  t = c(a(o.options, !0));
                return (
                  (t._handler = function (e) {
                    e.source === t.contentWindow &&
                      e.origin === Dropbox.baseUrl &&
                      p(e, null, o.options);
                  }),
                  Dropbox.addListener(window, "message", t._handler),
                  t
                );
              }),
              (Dropbox.cleanupWidget = function (e) {
                if (!e._handler) throw new Error("Invalid widget!");
                Dropbox.removeListener(window, "message", e._handler),
                  delete e._handler;
              });
          });
      },
      function (e, o, t) {
        const n =
          (this && this.__spreadArrays) ||
          function () {
            for (var e = 0, o = 0, t = arguments.length; o < t; o++)
              e += arguments[o].length;
            let n = Array(e),
              r = 0;
            for (o = 0; o < t; o++)
              for (let i = arguments[o], s = 0, a = i.length; s < a; s++, r++)
                n[r] = i[s];
            return n;
          };
        Object.defineProperty(o, "__esModule", { value: !0 });
        const r = window.location.protocol + "//" + window.location.host,
          i = (function () {
            function e(o) {
              if (
                ((this.origin = r),
                (this.sendMessage = function (e) {}),
                (this.state = {}),
                (this.options = o),
                !this.options)
              )
                throw new Error("options must be provided");
              if (!this.options.appKey)
                throw new Error("appKey must be provided");
              e.validateOnError(this.options.onError);
            }
            return (
              (e.validateOnError = function (e) {
                if (e && "function" != typeof e)
                  throw new Error("onError must be a function");
              }),
              (e.prototype.setOnError = function (o) {
                e.validateOnError(o), (this.options.onError = o);
              }),
              (e.prototype.hasOnCloseDialogMessage = function () {
                return void 0 !== this.onCloseDialogMessage;
              }),
              (e.prototype.setOnCloseDialogMessage = function (e) {
                if ("function" != typeof e)
                  throw new Error("onCloseDialogMessage must be a function");
                this.onCloseDialogMessage = e;
              }),
              (e.prototype.sendState = function () {
                this.sendMessage({ method: "state", params: this.state });
              }),
              (e.prototype.url = function () {
                const e = n(
                  [
                    { key: "app_key", value: this.options.appKey },
                    { key: "origin", value: this.origin },
                  ],
                  this.urlParams()
                )
                  .map(function (e) {
                    return (
                      encodeURIComponent(e.key) +
                      "=" +
                      encodeURIComponent(e.value)
                    );
                  })
                  .join("&");
                return { pathname: this.urlPathname(), search: "?" + e };
              }),
              (e.prototype.windowDimensions = function () {
                return { width: 735, height: 552 };
              }),
              (e.prototype.handleMessage = function (e) {
                switch (
                  (void 0 !== e.sequenceNumber &&
                    this.sendMessage({
                      method: "ack",
                      sequenceNumber: e.sequenceNumber,
                    }),
                  e.method)
                ) {
                  case "origin_request":
                    this.sendMessage({ method: "origin" });
                    break;
                  case "ready":
                    this.sendState();
                    break;
                  case "error":
                    this.options.onError && this.options.onError(e.params);
                    break;
                  case "close_dialog":
                    this.onCloseDialogMessage && this.onCloseDialogMessage(),
                      (this.onCloseDialogMessage = void 0);
                }
              }),
              e
            );
          })();
        o.Dropin = i;
      },
      function (e, o, t) {
        Object.defineProperty(o, "__esModule", { value: !0 });
        const n = t(0),
          r = ["text", "documents", "images", "video", "audio"];
        function i(e, o) {
          null != o
            ? (o.innerHTML = "")
            : ((o = document.createElement("a")).href = "#"),
            (o.className += " dropbox-dropin-btn"),
            Dropbox.isBrowserSupported()
              ? (o.className += " dropbox-dropin-default")
              : (o.className += " dropbox-dropin-disabled");
          const t = document.createElement("span");
          return (
            (t.className = "dropin-btn-status"),
            o.appendChild(t),
            (e = document.createTextNode(e)),
            o.appendChild(e),
            o
          );
        }
        function s(e) {
          return e.replace(/\/+$/g, "").split("/").pop();
        }
        function a(e) {
          const o = document.createElement("a");
          return (o.href = e), s(o.pathname);
        }
        (o.genericDropins = { init: function () {} }),
          (o.createDropinButton = i),
          (o.filenameFromPath = s),
          (o.initModule = function () {
            let e;
            n.initModule(),
              null == Dropbox.appKey &&
                (Dropbox.appKey =
                  null === (e = document.getElementById("dropboxjs")) ||
                  void 0 === e
                    ? void 0
                    : e.getAttribute("data-app-key")),
              (Dropbox.init = function (e) {
                null != e.appKey && (Dropbox.appKey = e.appKey);
              });
            const t = function (e) {
              let o, t, n;
              if ("string" == typeof e[0])
                (n = e.shift()),
                  (o = "string" == typeof e[0] ? e.shift() : a(n)),
                  ((t = e.shift() || {}).files = [{ url: n, filename: o }]);
              else {
                if (null == (t = e.shift()))
                  throw new Error("Missing arguments. See documentation.");
                if (
                  (null == t.files || !t.files.length) &&
                  "function" != typeof t.files
                )
                  throw new Error("Missing files. See documentation.");
                if (null != t.fetch_urls_fn) {
                  if ("function" != typeof t.fetch_urls_fn)
                    throw new Error(
                      "fetch_urls_fn must be a function if supplied.  See documentation."
                    );
                  t._fetch_url_on_save = !0;
                }
                for (let r = 0; r < t.files.length; r++) {
                  const i = t.files[r];
                  if (
                    "function" == typeof i.url &&
                    ((t._fetch_url_on_save = !0),
                    (t.fetch_urls_fn = i.url),
                    (i.url = null),
                    r > 0)
                  )
                    throw new Error(
                      "Old style url as callback is only supported for single files."
                    );
                  i.filename || (i.filename = a(i.url));
                }
              }
              return t;
            };
            Dropbox.save = function () {
              for (var e = [], o = 0; o < arguments.length; o++)
                e[o] = arguments[o];
              const r = t(e);
              if (Dropbox.isBrowserSupported()) {
                if (
                  ((r._popup = !0),
                  "object" != typeof r.files || !r.files.length)
                )
                  throw new Error(
                    "The object passed in must have a 'files' property that contains a list of objects. See documentation."
                  );
                if (r.iframe && !r.windowName)
                  throw new Error(
                    "Dropbox.save does not yet support creating its own iframe.                       windowName must be provided when the iframe option is present."
                  );
                for (let i = 0, a = r.files; i < a.length; i++) {
                  const l = a[i];
                  if (r._fetch_url_on_save) {
                    if (r.fetch_urls_fn) {
                      if (null != l.url)
                        throw new Error(
                          "You passed in a 'fetch_urls_fn' option to specify the file URLs.  Don't include individual URLs in each file objects."
                        );
                    } else if ("function" != typeof l.url)
                      throw new Error(
                        "File urls should be all urls, or a single file with function. See documentation."
                      );
                  } else if ("string" != typeof l.url)
                    throw new Error(
                      "File urls to download incorrectly configured. Each file must have a url. See documentation."
                    );
                }
                const c = n.popupDimensionsString(735, 670);
                return s(n.saverUrl(r.iframe), c, r).window;
              }
              alert("Your browser does not support the Dropbox Saver");
            };
            var s = function (e, o, t) {
                var r = function () {
                    a.closed ||
                      (a.close(),
                      a.postMessage(
                        JSON.stringify({ method: "close" }),
                        Dropbox.baseUrl
                      )),
                      Dropbox.removeListener(window, "message", i),
                      clearInterval(l);
                  },
                  i = function (e) {
                    (e.source !== a &&
                      e.source !==
                        (void 0 !== n.ieframe && null !== n.ieframe
                          ? n.ieframe.contentWindow
                          : void 0)) ||
                      n.handleMessageEvent(e, r, t);
                  },
                  s = t.iframe ? "" : o + ",resizable,scrollbars",
                  a = window.open(e, t.windowName || "dropbox", s);
                if (!a)
                  throw new Error(
                    "Failed to open/load the window. Dropbox.choose and Dropbox.save should only be called from within a user-triggered event handler such as a tap or click event."
                  );
                a.focus();
                var l = setInterval(function () {
                  (function () {
                    try {
                      return a.closed;
                    } catch (e) {}
                  })() &&
                    (r(),
                    "function" == typeof t.cancel &&
                      t.cancel(n.last_navigation));
                }, 100);
                return (
                  Dropbox.addListener(window, "message", i),
                  { window: a, onClose: r }
                );
              },
              l = function (e) {
                null == e.success &&
                  n.__guardMethod__(console, "warn", function (e) {
                    return e.warn(
                      "You must provide a success callback to the Chooser to see the files that the user selects"
                    );
                  }),
                  void 0 === e.fileselect ||
                    Boolean(e.fileselect) ||
                    Boolean(e.folderselect) ||
                    n.__guardMethod__(console, "error", function (e) {
                      return e.error(
                        "You must enable either fileselect or folderselect on the Chooser so the user can select something"
                      );
                    });
                const o = function () {
                  return (
                    n.__guardMethod__(console, "warn", function (e) {
                      return e.warn(
                        "The provided list of extensions or file types is not valid. See Chooser documentation: https://www.dropbox.com/developers/dropins/chooser/js"
                      );
                    }),
                    n.__guardMethod__(console, "warn", function (e) {
                      return e.warn(
                        "Available file types are: " + r.join(", ")
                      );
                    }),
                    delete e.extensions
                  );
                };
                if (null != e.extensions && null != Array.isArray)
                  if (Array.isArray(e.extensions))
                    for (let t = 0, i = e.extensions; t < i.length; t++) {
                      const s = i[t];
                      s.match(/^\.[\.\w$#&+@!()\-'`_~]+$/) ||
                        -1 !== r.indexOf(s) ||
                        o();
                    }
                  else o();
                return (
                  void 0 !== e.sizeLimit &&
                    "number" != typeof e.sizeLimit &&
                    e.sizeLimit <= 0 &&
                    n.__guardMethod__(console, "error", function (e) {
                      return e.error(
                        "The sizeLimit option, if provided, must be a positive number"
                      );
                    }),
                  e
                );
              },
              c = function (e) {
                if (Dropbox.isBrowserSupported()) {
                  let o,
                    t,
                    r = n.createChooserSession(e);
                  if (e.iframe && !e.windowName) {
                    const i =
                      ((o = n.chooserUrl(e, !0)),
                      ((t = document.createElement("iframe")).src = o),
                      (t.style.display = "block"),
                      (t.style.backgroundColor = "white"),
                      (t.style.border = "none"),
                      t);
                    (i.style.width = "735px"),
                      (i.style.height = "552px"),
                      (i.style.margin = "125px auto 0 auto"),
                      (i.style.border = "1px solid #ACACAC"),
                      (i.style.boxShadow = "rgba(0, 0, 0, .2) 0px 4px 16px");
                    const a = document.createElement("div");
                    (a.style.position = "fixed"),
                      (a.style.left = a.style.right = a.style.top = a.style.bottom =
                        "0"),
                      (a.style.zIndex = "1000"),
                      (a.style.backgroundColor = "rgba(160, 160, 160, 0.2)"),
                      a.appendChild(i),
                      document.body.appendChild(a);
                    var l = function (e) {
                      e.source === i.contentWindow &&
                        ((r.onClose = function () {
                          document.body.removeChild(a),
                            Dropbox.removeListener(window, "message", l);
                        }),
                        n.handleMessageEvent(e, r.onClose, r.options));
                    };
                    Dropbox.addListener(window, "message", l);
                  } else {
                    const c = n.popupDimensionsString(735, 552);
                    r.onClose = s(
                      n.chooserUrl(r.options, r.options.iframe),
                      c,
                      r.options
                    ).onClose;
                  }
                } else
                  alert("Your browser does not support the Dropbox Chooser");
              };
            (Dropbox.choose = function (e) {
              null == e && (e = {}), (e = l(e)), c(e);
            }),
              (Dropbox.cancelChooser = function () {
                n.currentChooserSession &&
                  (n.currentChooserSession.onClose &&
                    n.currentChooserSession.onClose(),
                  n.currentChooserSession.options.cancel &&
                    n.currentChooserSession.options.cancel(n.last_navigation));
              }),
              (Dropbox.isBrowserSupported = function () {
                const e = (function () {
                  for (
                    let e = 0, o = [/IEMobile\/(7|8|9|10)\./, /BB10;/, /CriOS/];
                    e < o.length;
                    e++
                  )
                    if (o[e].test(navigator.userAgent)) return !1;
                  return (
                    "undefined" != typeof JSON &&
                    null !== JSON &&
                    null != window.postMessage &&
                    null != window.addEventListener &&
                    !/MSIE [7-9]/.test(navigator.userAgent)
                  );
                })();
                return (
                  (Dropbox.isBrowserSupported = function () {
                    return e;
                  }),
                  e
                );
              }),
              (Dropbox.createChooseButton = function (e) {
                null == e && (e = {}), (e = l(e));
                const o = i("Choose from Dropbox");
                return (
                  Dropbox.addListener(o, "click", function (t) {
                    t.preventDefault(),
                      c({
                        success: function (t, n) {
                          (o.className =
                            "dropbox-dropin-btn dropbox-dropin-success"),
                            "function" == typeof e.success && e.success(t, n);
                        },
                        cancel: e.cancel,
                        linkType: e.linkType,
                        multiselect: e.multiselect,
                        fileselect: e.fileselect,
                        folderselect: e.folderselect,
                        extensions: e.extensions,
                        sizeLimit: e.sizeLimit,
                        iframe: e.iframe,
                        requiredPermissions: e.requiredPermissions,
                        _trigger: "button",
                      });
                  }),
                  o
                );
              }),
              (Dropbox.createSaveButton = function () {
                for (var e = [], o = 0; o < arguments.length; o++)
                  e[o] = arguments[o];
                let n = t(e),
                  r = e.shift();
                return (
                  (r = i("Save to Dropbox", r)),
                  Dropbox.addListener(r, "click", function (e) {
                    if (
                      (e.preventDefault(),
                      r.className.indexOf("dropbox-dropin-error") >= 0 ||
                        r.className.indexOf("dropbox-dropin-default") >= 0 ||
                        r.className.indexOf("dropbox-dropin-disabled") >= 0)
                    ) {
                      const o =
                        ("function" == typeof n.files ? n.files() : void 0) ||
                        n.files;
                      if (!(null != o ? o.length : void 0))
                        return (
                          (r.className =
                            "dropbox-dropin-btn dropbox-dropin-error"),
                          void (
                            "function" == typeof n.error &&
                            n.error("Missing files")
                          )
                        );
                      Dropbox.save({
                        files: o,
                        success: function () {
                          (r.className =
                            "dropbox-dropin-btn dropbox-dropin-success"),
                            "function" == typeof n.success && n.success();
                        },
                        progress: function (e) {
                          (r.className =
                            "dropbox-dropin-btn dropbox-dropin-progress"),
                            "function" == typeof n.progress && n.progress(e);
                        },
                        cancel: function () {
                          "function" == typeof n.cancel && n.cancel();
                        },
                        error: function (e) {
                          (r.className =
                            "dropbox-dropin-btn dropbox-dropin-error"),
                            "function" == typeof n.error && n.error(e);
                        },
                      });
                    }
                  }),
                  r
                );
              });
            const u = function (e, o) {
                return (
                  "  background: " +
                  e +
                  ";\n  background: -moz-linear-gradient(top, " +
                  e +
                  " 0%, " +
                  o +
                  " 100%);\n  background: -webkit-linear-gradient(top, " +
                  e +
                  " 0%, " +
                  o +
                  " 100%);\n  background: linear-gradient(to bottom, " +
                  e +
                  " 0%, " +
                  o +
                  " 100%);\n  filter: progid:DXImageTransform.Microsoft.gradient(startColorstr='" +
                  e +
                  "', endColorstr='" +
                  o +
                  "',GradientType=0);  "
                );
              },
              p = document.createElement("style");
            p.type = "text/css";
            const d =
              '  @-webkit-keyframes rotate {\n    from  { -webkit-transform: rotate(0deg); }\n    to   { -webkit-transform: rotate(360deg); }\n  }\n\n  @keyframes rotate {\n    from  { transform: rotate(0deg); }\n    to   { transform: rotate(360deg); }\n  }\n\n    .dropbox-dropin-btn, .dropbox-dropin-btn:link, .dropbox-dropin-btn:hover {\n      display: inline-block;\n      height: 14px;\n      font-family: "Lucida Grande", "Segoe UI", "Tahoma", "Helvetica Neue", "Helvetica", sans-serif;\n      font-size: 11px;\n      font-weight: 600;\n      color: #636363;\n      text-decoration: none;\n      padding: 1px 7px 5px 3px;\n      border: 1px solid #ebebeb;\n      border-radius: 2px;\n      border-bottom-color: #d4d4d4;\n      ' +
              u("#fcfcfc", "#f5f5f5") +
              "\n    }\n\n    .dropbox-dropin-default:hover, .dropbox-dropin-error:hover {\n      border-color: #dedede;\n      border-bottom-color: #cacaca;\n      " +
              u("#fdfdfd", "#f5f5f5") +
              "\n    }\n\n    .dropbox-dropin-default:active, .dropbox-dropin-error:active {\n      border-color: #d1d1d1;\n      box-shadow: inset 0 1px 1px rgba(0,0,0,0.1);\n    }\n\n    .dropbox-dropin-btn .dropin-btn-status {\n      display: inline-block;\n      width: 15px;\n      height: 14px;\n      vertical-align: bottom;\n      margin: 0 5px 0 2px;\n      background: transparent url('" +
              Dropbox.baseUrl +
              "/static/images/widgets/dbx-saver-status.png') no-repeat;\n      position: relative;\n      top: 2px;\n    }\n\n    .dropbox-dropin-default .dropin-btn-status {\n      background-position: 0px 0px;\n    }\n\n    .dropbox-dropin-progress .dropin-btn-status {\n      width: 18px;\n      margin: 0 4px 0 0;\n      background: url('" +
              Dropbox.baseUrl +
              "/static/images/widgets/dbx-progress.png') no-repeat center center;\n        -webkit-animation-name: rotate;\n        -webkit-animation-duration: 1.7s;\n        -webkit-animation-iteration-count: infinite;\n        -webkit-animation-timing-function: linear;\n      animation-name: rotate;\n      animation-duration: 1.7s;\n      animation-iteration-count: infinite;\n      animation-timing-function: linear;\n    }\n\n    .dropbox-dropin-success .dropin-btn-status {\n      background-position: -15px 0px;\n    }\n\n    .dropbox-dropin-disabled {\n      background: #e0e0e0;\n      border: 1px #dadada solid;\n      border-bottom: 1px solid #ccc;\n      box-shadow: none;\n    }\n\n    .dropbox-dropin-disabled .dropin-btn-status {\n      background-position: -30px 0px;\n    }\n\n    .dropbox-dropin-error .dropin-btn-status {\n      background-position: -45px 0px;\n    }\n\n  @media only screen and (-webkit-min-device-pixel-ratio: 1.4) {\n      .dropbox-dropin-btn .dropin-btn-status {\n        background-image: url('" +
              Dropbox.baseUrl +
              "/static/images/widgets/dbx-saver-status-2x.png');\n        background-size: 60px 14px;\n          -webkit-background-size: 60px 14px;\n      }\n\n      .dropbox-dropin-progress .dropin-btn-status {\n        background: url('" +
              Dropbox.baseUrl +
              "/static/images/widgets/dbx-progress-2x.png') no-repeat center center;\n        background-size: 20px 20px;\n          -webkit-background-size: 20px 20px;\n      }\n  }\n\n    .dropbox-saver:hover, .dropbox-chooser:hover {\n      text-decoration: none;\n      cursor: pointer;\n    }\n\n    .dropbox-chooser, .dropbox-dropin-btn {\n      line-height: 11px !important;\n      text-decoration: none !important;\n      box-sizing: content-box !important;\n        -webkit-box-sizing: content-box !important;\n        -moz-box-sizing: content-box !important;\n    }\n    ";
            p.styleSheet ? (p.styleSheet.cssText = d) : (p.textContent = d),
              document.getElementsByTagName("head")[0].appendChild(p),
              setTimeout(n.createIEFrame, 0);
            var f = function () {
              document.removeEventListener
                ? document.removeEventListener("DOMContentLoaded", f, !1)
                : document.detachEvent &&
                  document.detachEvent("onreadystatechange", f),
                n.createIEFrame(),
                o.genericDropins.init();
            };
            "interactive" === document.readyState ||
            "complete" === document.readyState
              ? setTimeout(f, 0)
              : document.addEventListener
              ? document.addEventListener("DOMContentLoaded", f, !1)
              : document.attachEvent("onreadystatechange", f);
          });
      },
      function (e, o, t) {
        let n,
          r =
            (this && this.__extends) ||
            ((n = function (e, o) {
              return (n =
                Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array &&
                  function (e, o) {
                    e.__proto__ = o;
                  }) ||
                function (e, o) {
                  for (const t in o) o.hasOwnProperty(t) && (e[t] = o[t]);
                })(e, o);
            }),
            function (e, o) {
              function t() {
                this.constructor = e;
              }
              n(e, o),
                (e.prototype =
                  null === o
                    ? Object.create(o)
                    : ((t.prototype = o.prototype), new t()));
            });
        Object.defineProperty(o, "__esModule", { value: !0 });
        const i = "https://www.dropbox.com/developers/dropins/chooser/js",
          s = ["text", "documents", "images", "video", "audio"],
          a = (function (e) {
            function o(o) {
              const t = e.call(this, o) || this;
              return (
                (t.validateOptions = function () {
                  const e = function (e, o) {
                    if (void 0 !== t.options[e] && typeof t.options[e] !== o)
                      throw new Error(
                        "The " + e + " option, if provided, must have type " + o
                      );
                  };
                  if (
                    (e("linkType", "string"),
                    e("_trigger", "string"),
                    void 0 !== t.options.extensions)
                  ) {
                    if (!(t.options.extensions instanceof Array))
                      throw new Error(
                        "The extensions option, if provided, must be an array"
                      );
                    for (
                      let o = 0, n = t.options.extensions;
                      o < n.length;
                      o++
                    ) {
                      const r = n[o];
                      if (
                        "string" != typeof r ||
                        (!r.match(/^\.[\.\w$#&+@!()\-'`_~]+$/) &&
                          -1 === s.indexOf(r))
                      )
                        throw new Error(
                          "The provided list of extensions or file types is not valid. See Chooser documentation: " +
                            i +
                            ". Available file types are: " +
                            s.join(", ")
                        );
                    }
                  }
                  if (
                    (e("multiselect", "boolean"),
                    e("iframe", "boolean"),
                    e("folderselect", "boolean"),
                    e("fileselect", "boolean"),
                    void 0 !== t.options.fileselect &&
                      !t.options.fileselect &&
                      !t.options.folderselect)
                  )
                    throw new Error(
                      "You must enable either fileselect or folderselect on the Chooser so the user can select something"
                    );
                  if (
                    (e("sizeLimit", "number"),
                    void 0 !== t.options.sizeLimit && t.options.sizeLimit <= 0)
                  )
                    throw new Error(
                      "The sizeLimit option, if provided, must be a positive number"
                    );
                  const a = t.options.initialNavigation;
                  if (void 0 !== a) {
                    if (void 0 !== a.mode && "string" != typeof a.mode)
                      throw new Error(
                        "The initialNavigation.mode option, if provided, must be a string"
                      );
                    if (void 0 !== a.role && "string" != typeof a.role)
                      throw new Error(
                        "The initialNavigation.role option, if provided, must be a string"
                      );
                    if (void 0 !== a.cursor && "string" != typeof a.cursor)
                      throw new Error(
                        "The initialNavigation.cursor option, if provided, must be a string"
                      );
                  }
                  if (
                    (e("initialViewType", "string"),
                    void 0 !== t.options.fields)
                  ) {
                    if (!(t.options.fields instanceof Array))
                      throw new Error(
                        "The fields option, if provided, must be an array"
                      );
                    for (let l = 0, c = t.options.fields; l < c.length; l++)
                      if ("string" != typeof c[l])
                        throw new Error(
                          "The fields option, if provided, must be an array of strings"
                        );
                  }
                  if (
                    (e("showSignOut", "boolean"),
                    e("version", "string"),
                    e("cl", "string"),
                    "function" != typeof t.options.onSuccess)
                  )
                    throw new Error(
                      "You must provide a success callback to the Chooser to see the files that the user selects"
                    );
                  e("onReady", "function"),
                    e("onCancel", "function"),
                    e("onError", "function"),
                    e("onResize", "function");
                }),
                t.validateOptions(),
                t
              );
            }
            return (
              r(o, e),
              (o.prototype.urlParams = function () {
                const e = [],
                  o = function (o, t) {
                    void 0 !== t && e.push({ key: o, value: "" + t });
                  };
                return (
                  o("link_type", this.options.linkType),
                  e.push({
                    key: "trigger",
                    value: this.options._trigger || "js",
                  }),
                  void 0 !== this.options.extensions &&
                    e.push({
                      key: "extensions",
                      value: this.options.extensions.join(" "),
                    }),
                  o("multiselect", this.options.multiselect),
                  o("iframe", this.options.iframe),
                  o("folderselect", this.options.folderselect),
                  o("fileselect", this.options.fileselect),
                  o("size_limit", this.options.sizeLimit),
                  void 0 !== this.options.initialNavigation &&
                    (o(
                      "initial_navigation_mode",
                      this.options.initialNavigation.mode
                    ),
                    o(
                      "initial_navigation_role",
                      this.options.initialNavigation.role
                    ),
                    o(
                      "initial_navigation_cursor",
                      this.options.initialNavigation.cursor
                    )),
                  o("initial_view_type", this.options.initialViewType),
                  void 0 !== this.options.fields &&
                    e.push({
                      key: "fields",
                      value: this.options.fields.join(" "),
                    }),
                  o("show_sign_out", this.options.showSignOut),
                  o("version", this.options.version),
                  o("cl", this.options.cl),
                  e
                );
              }),
              (o.prototype.urlPathname = function () {
                return "/chooser";
              }),
              (o.prototype.close = function () {
                "function" == typeof this.onCloseDialogMessage &&
                  this.onCloseDialogMessage();
              }),
              (o.prototype.handleMessage = function (o) {
                switch ((e.prototype.handleMessage.call(this, o), o.method)) {
                  case "ready":
                    void 0 !== this.options.onReady && this.options.onReady();
                    break;
                  case "files_selected":
                    this.close(),
                      this.options.onSuccess(o.params, this.lastNavigation);
                    break;
                  case "cursor_changed":
                    this.lastNavigation = { cursor: o.params };
                    break;
                  case "close_dialog":
                    void 0 !== this.options.onCancel &&
                      this.options.onCancel(this.lastNavigation);
                    break;
                  case "resize":
                    void 0 !== this.options.onResize &&
                      this.options.onResize(o.params);
                    break;
                  case "error":
                    this.close();
                    break;
                  case "_debug_log":
                    void 0 !== console &&
                      null !== console &&
                      console.log(o.params.msg);
                }
              }),
              o
            );
          })(t(1).Dropin);
        o.BaseChooser = a;
        const l = (function (e) {
          function o(o) {
            return e.call(this, o) || this;
          }
          return r(o, e), o;
        })(a);
        o.Chooser = l;
      },
      ,
      ,
      function (e, o, t) {
        Object.defineProperty(o, "__esModule", { value: !0 }),
          t(7).initModule(),
          (o.Dropbox = window.Dropbox);
      },
      function (e, o, t) {
        var n =
          (this && this.__assign) ||
          function () {
            return (n =
              Object.assign ||
              function (e) {
                for (var o, t = 1, n = arguments.length; t < n; t++)
                  for (const r in (o = arguments[t]))
                    Object.prototype.hasOwnProperty.call(o, r) && (e[r] = o[r]);
                return e;
              }).apply(this, arguments);
          };
        Object.defineProperty(o, "__esModule", { value: !0 });
        const r = t(8),
          i = t(2),
          s = t(3),
          a = t(10),
          l = t(11),
          c = t(12),
          u = t(13);
        o.initModule = function () {
          i.initModule(), (Dropbox.VERSION = "2");
          const e = new r.BrowserEnvironment();
          (Dropbox.mount = e.mount.bind(e)),
            (Dropbox.openWindow = e.openWindow.bind(e));
          const o = e.remove.bind(e);
          (Dropbox.unmount = o),
            (Dropbox.closeWindow = o),
            (Dropbox.Mover = c.Mover),
            (Dropbox.Previewer = u.Previewer),
            (Dropbox.Chooser = s.Chooser),
            (Dropbox.ZoomChooser = a.ZoomChooser),
            (Dropbox.embed = function (e, o) {
              Dropbox.appKey && (e = n(n({}, e), { appKey: Dropbox.appKey }));
              const t = new l.Embed(e);
              return Dropbox.mount(t, o), t;
            }),
            (i.genericDropins.init = function () {
              for (
                let e = document.getElementsByTagName("a"), o = e.length - 1;
                o >= 0;
                o--
              ) {
                const t = e[o],
                  n = (t.getAttribute("class") || "").split(" ");
                n.indexOf("dropbox-saver") >= 0
                  ? (function (e) {
                      Dropbox.createSaveButton(
                        {
                          files: function () {
                            return [
                              {
                                url: e.getAttribute("data-url") || e.href,
                                filename:
                                  e.getAttribute("data-filename") ||
                                  i.filenameFromPath(e.pathname),
                              },
                            ];
                          },
                        },
                        e
                      );
                    })(t)
                  : n.indexOf("dropbox-embed") >= 0 &&
                    (function (e) {
                      const o = e.getAttribute("data-url") || e.href;
                      if (o && e.parentElement) {
                        const t = e.getAttribute("data-file-zoom") || void 0,
                          n = e.getAttribute("data-folder-view") || void 0,
                          r =
                            e.getAttribute("data-folder-header-size") || void 0,
                          i = document.createElement("div");
                        i.classList.add("dropbox-embed-container"),
                          (i.style.height =
                            e.getAttribute("data-height") || "100%"),
                          (i.style.width =
                            e.getAttribute("data-width") || "100%"),
                          e.parentElement.replaceChild(i, e),
                          Dropbox.embed(
                            {
                              link: o,
                              file: { zoom: t },
                              folder: { view: n, headerSize: r },
                            },
                            i
                          );
                      }
                    })(t);
              }
            });
        };
      },
      function (e, o, t) {
        Object.defineProperty(o, "__esModule", { value: !0 });
        const n = t(9);
        function r(e) {
          return function () {
            for (let o = 0, t = e; o < t.length; o++) (0, t[o])();
          };
        }
        o.TARGET_ORIGIN = "https://www.dropbox.com";
        const i = function () {},
          s = (function () {
            function e() {
              const e = this;
              (this.activeDropins = []),
                (this.deleteActiveDropin = function (o) {
                  return function () {
                    const t = e.activeDropins.indexOf(o);
                    -1 !== t && e.activeDropins.splice(t, 1);
                  };
                }),
                (this.openWindow = function (t) {
                  e.throwIfAlreadyActive(t);
                  let s = r([]);
                  try {
                    const a = t.url(),
                      l = "" + o.TARGET_ORIGIN + a.pathname + a.search,
                      c = n.PopupEnvironment.open(
                        l,
                        t.windowDimensions(),
                        function () {
                          s();
                        }
                      );
                    s = r([c.stopInterval, s]);
                    const u = e.attach(t, c.messagingWindow);
                    s = r([
                      function () {
                        t.sendMessage({ method: "close" });
                      },
                      (s = r([u, s])),
                    ]);
                    const p = { dropin: t, cleanup: i };
                    e.activeDropins.push(p),
                      (s = r([e.deleteActiveDropin(p), s])),
                      t.hasOnCloseDialogMessage() ||
                        t.setOnCloseDialogMessage(s),
                      (p.cleanup = s);
                  } catch (e) {
                    throw (s(), e);
                  }
                });
            }
            return (
              (e.prototype.mount = function (e, t) {
                if (!e) throw new Error("Dropbox component must be provided");
                if (!t) throw new Error("Container element must be provided");
                this.throwIfAlreadyActive(e);
                let n = r([]);
                try {
                  const s = this.createIframe();
                  n = r([
                    this.attach(e, function () {
                      if (!s.contentWindow)
                        throw new Error(
                          "iframe does not contain a contentWindow"
                        );
                      return s.contentWindow;
                    }),
                    n,
                  ]);
                  const a = e.url();
                  (s.src =
                    "" +
                    o.TARGET_ORIGIN +
                    a.pathname +
                    a.search +
                    "&iframe=true"),
                    (s.scrolling = "no"),
                    t.appendChild(s),
                    (n = r([
                      function () {
                        t.removeChild(s);
                      },
                      n,
                    ]));
                  const l = { dropin: e, cleanup: i };
                  this.activeDropins.push(l),
                    (n = r([this.deleteActiveDropin(l), n])),
                    e.hasOnCloseDialogMessage() || e.setOnCloseDialogMessage(n),
                    (l.cleanup = n);
                } catch (e) {
                  throw (n(), e);
                }
              }),
              (e.prototype.remove = function (e) {
                if (!e) throw new Error("Dropbox component must be provided");
                for (let o = 0, t = this.activeDropins; o < t.length; o++) {
                  const n = t[o];
                  if (n.dropin === e) {
                    n.cleanup();
                    break;
                  }
                }
              }),
              (e.prototype.throwIfAlreadyActive = function (e) {
                for (let o = 0, t = this.activeDropins; o < t.length; o++)
                  if (t[o].dropin === e)
                    throw new Error("Component is already in use");
              }),
              (e.prototype.attach = function (e, t) {
                e.sendMessage = function (e) {
                  t().postMessage(JSON.stringify(e), o.TARGET_ORIGIN);
                };
                const n = function (n) {
                  if (n.source === t() && n.origin === o.TARGET_ORIGIN) {
                    let r;
                    try {
                      r = JSON.parse(n.data);
                    } catch (e) {
                      return;
                    }
                    "object" == typeof r &&
                      "string" == typeof r.method &&
                      e.handleMessage(r);
                  }
                };
                return (
                  window.addEventListener("message", n),
                  function () {
                    (e.sendMessage = i),
                      window.removeEventListener("message", n);
                  }
                );
              }),
              (e.prototype.createIframe = function () {
                const e = window.document.createElement("iframe");
                return (
                  (e.style.height = "100%"),
                  (e.style.width = "100%"),
                  (e.style.border = "none"),
                  (e.allowFullscreen = !0),
                  e
                );
              }),
              e
            );
          })();
        o.BrowserEnvironment = s;
      },
      function (e, o, t) {
        Object.defineProperty(o, "__esModule", { value: !0 });
        const n = t(0),
          r = (function () {
            function e(e, o) {
              const t = this;
              (this.popupWindow = e),
                (this.onClose = o),
                (this.stopInterval = function () {
                  clearInterval(t.intervalId);
                }),
                (this.isWindowClosedByUser = function () {
                  try {
                    return t.popupWindow.closed;
                  } catch (e) {}
                  return !1;
                }),
                (this.messagingWindow = function () {
                  return void 0 !== n.ieframe && null !== n.ieframe
                    ? n.ieframe.contentWindow
                    : t.popupWindow;
                }),
                (this.handleInterval = function () {
                  t.isWindowClosedByUser() && (t.onClose(), t.stopInterval());
                }),
                (this.intervalId = setInterval(this.handleInterval, 100));
            }
            return (
              (e.open = function (o, t, r) {
                const i =
                    n.popupDimensionsString(t.width, t.height) +
                    ",resizable,scrollbars",
                  s = window.open(o, "_blank", i);
                if (null === s)
                  throw new Error(
                    "Failed to open the window. Dropbox popups may only be attached to a user-triggered event handler such as a tap or click event."
                  );
                return s.focus(), new e(s, r);
              }),
              e
            );
          })();
        o.PopupEnvironment = r;
      },
      function (e, o, t) {
        let n,
          r =
            (this && this.__extends) ||
            ((n = function (e, o) {
              return (n =
                Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array &&
                  function (e, o) {
                    e.__proto__ = o;
                  }) ||
                function (e, o) {
                  for (const t in o) o.hasOwnProperty(t) && (e[t] = o[t]);
                })(e, o);
            }),
            function (e, o) {
              function t() {
                this.constructor = e;
              }
              n(e, o),
                (e.prototype =
                  null === o
                    ? Object.create(o)
                    : ((t.prototype = o.prototype), new t()));
            });
        Object.defineProperty(o, "__esModule", { value: !0 });
        const i = (function (e) {
          function o(o) {
            const t = e.call(this, o) || this;
            return t.validateNoOptionLinkType(), t;
          }
          return (
            r(o, e),
            (o.prototype.validateNoOptionLinkType = function () {
              if (void 0 !== this.options.linkType)
                throw new Error(
                  "The ZoomChooser SDK does not accept a linkType."
                );
            }),
            (o.prototype.urlParams = function () {
              const o = e.prototype.urlParams.call(this);
              return o.push({ key: "is_zoom_chooser", value: "true" }), o;
            }),
            o
          );
        })(t(3).BaseChooser);
        o.ZoomChooser = i;
      },
      function (e, o, t) {
        let n,
          r =
            (this && this.__extends) ||
            ((n = function (e, o) {
              return (n =
                Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array &&
                  function (e, o) {
                    e.__proto__ = o;
                  }) ||
                function (e, o) {
                  for (const t in o) o.hasOwnProperty(t) && (e[t] = o[t]);
                })(e, o);
            }),
            function (e, o) {
              function t() {
                this.constructor = e;
              }
              n(e, o),
                (e.prototype =
                  null === o
                    ? Object.create(o)
                    : ((t.prototype = o.prototype), new t()));
            }),
          i =
            (this && this.__spreadArrays) ||
            function () {
              for (var e = 0, o = 0, t = arguments.length; o < t; o++)
                e += arguments[o].length;
              let n = Array(e),
                r = 0;
              for (o = 0; o < t; o++)
                for (let i = arguments[o], s = 0, a = i.length; s < a; s++, r++)
                  n[r] = i[s];
              return n;
            };
        Object.defineProperty(o, "__esModule", { value: !0 });
        const s = (function (e) {
          function o(o) {
            let t, n, r;
            void 0 === o && (o = {});
            const i = e.call(this, o) || this;
            return (
              (i.state = {
                link: i.options.link,
                file_zoom:
                  null === (t = i.options.file) || void 0 === t
                    ? void 0
                    : t.zoom,
                folder_view:
                  null === (n = i.options.folder) || void 0 === n
                    ? void 0
                    : n.view,
                folder_headerSize:
                  null === (r = i.options.folder) || void 0 === r
                    ? void 0
                    : r.headerSize,
              }),
              i
            );
          }
          return (
            r(o, e),
            (o.prototype.urlPathname = function () {
              return "/dropins/embed";
            }),
            (o.prototype.urlParams = function () {
              return i(
                this.state.link
                  ? [{ key: "link", value: this.state.link }]
                  : [],
                this.state.file_zoom
                  ? [{ key: "file_zoom", value: this.state.file_zoom }]
                  : [],
                this.state.folder_view
                  ? [{ key: "folder_view", value: this.state.folder_view }]
                  : [],
                this.state.folder_headerSize
                  ? [
                      {
                        key: "folder_header_size",
                        value: this.state.folder_headerSize,
                      },
                    ]
                  : []
              );
            }),
            (o.prototype.handleMessage = function (o) {
              switch (o.method) {
                case "update_size":
                  this.options.onSizeChanged &&
                    this.options.onSizeChanged(o.params);
                  break;
                default:
                  e.prototype.handleMessage.call(this, o);
              }
            }),
            o
          );
        })(t(1).Dropin);
        o.Embed = s;
      },
      function (e, o, t) {
        var n,
          r =
            (this && this.__extends) ||
            ((n = function (e, o) {
              return (n =
                Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array &&
                  function (e, o) {
                    e.__proto__ = o;
                  }) ||
                function (e, o) {
                  for (const t in o) o.hasOwnProperty(t) && (e[t] = o[t]);
                })(e, o);
            }),
            function (e, o) {
              function t() {
                this.constructor = e;
              }
              n(e, o),
                (e.prototype =
                  null === o
                    ? Object.create(o)
                    : ((t.prototype = o.prototype), new t()));
            }),
          i =
            (this && this.__assign) ||
            function () {
              return (i =
                Object.assign ||
                function (e) {
                  for (var o, t = 1, n = arguments.length; t < n; t++)
                    for (const r in (o = arguments[t]))
                      Object.prototype.hasOwnProperty.call(o, r) &&
                        (e[r] = o[r]);
                  return e;
                }).apply(this, arguments);
            };
        Object.defineProperty(o, "__esModule", { value: !0 });
        const s = (function (e) {
          function o(t) {
            const n = e.call(this, t) || this;
            return (
              o.validateOnSuccess(n.options.onSuccess),
              o.validateOnCancel(n.options.onCancel),
              (n.state = {
                entries: n.options.entries,
                initialFolderSelection: n.options.initialFolderSelection,
              }),
              (n.onSuccess = n.options.onSuccess),
              (n.onCancel = n.options.onCancel),
              (n.onFolderSelected = n.options.onFolderSelected),
              n
            );
          }
          return (
            r(o, e),
            (o.validateOnSuccess = function (e) {
              if (e && "function" != typeof e)
                throw new Error("onSuccess must be a function");
            }),
            (o.validateOnCancel = function (e) {
              if (e && "function" != typeof e)
                throw new Error("onCancel must be a function");
            }),
            (o.validateOnFolderSelected = function (e) {
              if (e && "function" != typeof e)
                throw new Error("onFolderSelected must be a function");
            }),
            (o.prototype.urlParams = function () {
              return [
                { key: "account_id", value: this.options.accountId },
                {
                  key: "initial_folder_selection",
                  value: this.state.initialFolderSelection,
                },
                {
                  key: "show_folders_only",
                  value: this.options.showFoldersOnly
                    ? String(this.options.showFoldersOnly)
                    : "",
                },
              ];
            }),
            (o.prototype.urlPathname = function () {
              return "/dropins/mover";
            }),
            (o.prototype.handleMessage = function (o) {
              switch ((e.prototype.handleMessage.call(this, o), o.method)) {
                case "success":
                  this.onSuccess && this.onSuccess();
                  break;
                case "cancel":
                  this.onCancel && this.onCancel();
                  break;
                case "folder_selected":
                  this.onFolderSelected && this.onFolderSelected(o.params.path);
              }
            }),
            (o.prototype.setEntries = function (e) {
              (this.state = i(i({}, this.state), { entries: e })),
                this.sendState();
            }),
            (o.prototype.setInitialFolderSelection = function (e) {
              (this.state = i(i({}, this.state), {
                initialFolderSelection: e,
              })),
                this.sendState();
            }),
            (o.prototype.setOnSuccess = function (e) {
              o.validateOnSuccess(e), (this.onSuccess = e);
            }),
            (o.prototype.setOnCancel = function (e) {
              o.validateOnCancel(e), (this.onCancel = e);
            }),
            (o.prototype.setOnFolderSelected = function (e) {
              o.validateOnFolderSelected(e), (this.onFolderSelected = e);
            }),
            o
          );
        })(t(1).Dropin);
        o.Mover = s;
      },
      function (e, o, t) {
        var n,
          r =
            (this && this.__extends) ||
            ((n = function (e, o) {
              return (n =
                Object.setPrototypeOf ||
                ({ __proto__: [] } instanceof Array &&
                  function (e, o) {
                    e.__proto__ = o;
                  }) ||
                function (e, o) {
                  for (const t in o) o.hasOwnProperty(t) && (e[t] = o[t]);
                })(e, o);
            }),
            function (e, o) {
              function t() {
                this.constructor = e;
              }
              n(e, o),
                (e.prototype =
                  null === o
                    ? Object.create(o)
                    : ((t.prototype = o.prototype), new t()));
            }),
          i =
            (this && this.__assign) ||
            function () {
              return (i =
                Object.assign ||
                function (e) {
                  for (var o, t = 1, n = arguments.length; t < n; t++)
                    for (const r in (o = arguments[t]))
                      Object.prototype.hasOwnProperty.call(o, r) &&
                        (e[r] = o[r]);
                  return e;
                }).apply(this, arguments);
            },
          s =
            (this && this.__spreadArrays) ||
            function () {
              for (var e = 0, o = 0, t = arguments.length; o < t; o++)
                e += arguments[o].length;
              let n = Array(e),
                r = 0;
              for (o = 0; o < t; o++)
                for (let i = arguments[o], s = 0, a = i.length; s < a; s++, r++)
                  n[r] = i[s];
              return n;
            };
        Object.defineProperty(o, "__esModule", { value: !0 });
        const a = (function (e) {
          function o(o) {
            void 0 === o && (o = {});
            const t = e.call(this, o) || this;
            return (
              (t.state = {
                accountId: t.options.accountId,
                hideAccount: t.options.hideAccount,
                link: t.options.link,
                cl: t.options.cl,
                view: t.options.initialView,
                openLinksWithSDK: !!t.options.onOpenLink,
              }),
              t
            );
          }
          return (
            r(o, e),
            (o.prototype.urlPathname = function () {
              return "/dropins/previewer";
            }),
            (o.prototype.urlParams = function () {
              return s(
                this.state.accountId
                  ? [{ key: "account_id", value: this.state.accountId }]
                  : [],
                this.state.hideAccount
                  ? [{ key: "hide_account", value: "true" }]
                  : [],
                this.state.view
                  ? [{ key: "initial_view", value: this.state.view }]
                  : [],
                this.state.cl ? [{ key: "cl", value: this.state.cl }] : [],
                this.state.link ? [{ key: "link", value: this.state.link }] : []
              );
            }),
            (o.prototype.handleMessage = function (o) {
              switch (o.method) {
                case "view_change":
                  var t = o.params.view;
                  (this.state = i(i({}, this.state), { view: t })),
                    this.options.onViewChange && this.options.onViewChange(t);
                  break;
                case "open_link":
                  this.options.onOpenLink &&
                    this.options.onOpenLink(o.params.link);
                  break;
                case "open_login":
                  this.options.onOpenLogin && this.options.onOpenLogin();
                  break;
                case "update_size":
                  this.options.onSizeChanged &&
                    this.options.onSizeChanged(o.params);
                  break;
                default:
                  e.prototype.handleMessage.call(this, o);
              }
            }),
            (o.prototype.setAccountId = function (e) {
              (this.state = i(i({}, this.state), { accountId: e })),
                this.sendState();
            }),
            (o.prototype.setHideAccount = function (e) {
              (this.state = i(i({}, this.state), { hideAccount: e })),
                this.sendState();
            }),
            (o.prototype.setLink = function (e) {
              (this.state = i(i({}, this.state), { link: e })),
                this.sendState();
            }),
            (o.prototype.setView = function (e) {
              (this.state = i(i({}, this.state), { view: e })),
                this.sendState();
            }),
            o
          );
        })(t(1).Dropin);
        o.Previewer = a;
      },
    ])
  );
}
/* eslint-enable */
