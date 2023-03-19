const configs = require(__dirname + "/src/configs/configs.js");
const nodeEnv = require(__dirname + "/src/common/constants/nodeEnv.js");
const path = require("path");
const inliner = require("@vgrid/sass-inline-svg");

const cloneDeep = require("lodash/cloneDeep");
const { createSecureHeaders } = require("next-secure-headers");

const isProdBuild = process.env.NODE_ENV === nodeEnv.PRODUCTION;

const PLAUSIBLE_URL = "https://plausible.io";

const TWITTER_URL = "https://cdn.syndication.twimg.com platform.twitter.com";

const WISTIA_URL = "https://fast.wistia.net/";

const SCRIPT_SRC = [
  "'self'",
  "'wasm-unsafe-eval'",
  PLAUSIBLE_URL,
  TWITTER_URL,
  WISTIA_URL,
];

const defaultSecureHeaders = {
  contentSecurityPolicy: {
    directives: {
      baseUri: "'self'",
      connectSrc: [
        "'self'",
        "sentry.prod.si.czi.technology",
        PLAUSIBLE_URL,
        configs.API_URL,
      ],
      defaultSrc: ["'self'"],
      fontSrc: ["'self'", "https://fonts.gstatic.com"],
      formAction: "'self'",
      frameAncestors: ["'none'"],
      frameSrc: [TWITTER_URL, WISTIA_URL],
      imgSrc: ["'self'", "data:"],
      manifestSrc: ["'self'"],
      mediaSrc: ["'self'"],
      objectSrc: ["'none'"],
      reportURI:
        configs.SENTRY_DEPLOYMENT_ENVIRONMENT &&
        "https://sentry.prod.si.czi.technology/api/167/security/?sentry_key=0432f3b3ceba4bc08d28dfb61fa29707&sentry_environment=" +
          configs.SENTRY_DEPLOYMENT_ENVIRONMENT,
      scriptSrc: isProdBuild ? SCRIPT_SRC : [...SCRIPT_SRC, "'unsafe-eval'"],
      styleSrc: ["'self'", "'unsafe-inline'", "https://fonts.googleapis.com"],
      upgradeInsecureRequests: true,
      workerSrc: true,
    },
  },
};

// unsafe-eval is required for next-mdx-remote
const docSiteScriptSrc = [...SCRIPT_SRC, "'unsafe-eval'"];
// Required for google slides iframe
const docSiteFrameSrc = ["https://docs.google.com"];
const docSiteSecureHeaders = cloneDeep(defaultSecureHeaders);
docSiteSecureHeaders.contentSecurityPolicy.directives.scriptSrc =
  docSiteScriptSrc;
docSiteSecureHeaders.contentSecurityPolicy.directives.frameSrc =
  docSiteFrameSrc;

module.exports = {
  compiler: { emotion: true },
  eslint: { dirs: ["doc-site", "pages", "components", "lib"] },
  async generateBuildId() {
    // Return null to allow next.js to fallback to default behavior
    // if COMMIT_SHA env is missing or empty.
    return process.env.COMMIT_SHA || null;
  },

  headers() {
    return [
      {
        headers: createSecureHeaders(defaultSecureHeaders),
        source: "/(.*)",
      },
      {
        headers: createSecureHeaders(docSiteSecureHeaders),
        source: `/docs/:slug*`,
      },
    ];
  },
  async redirects() {
    return [
      {
        destination: "/docs/01__CellxGene",
        permanent: true,
        source: "/docs",
      },
      {
        destination: "/",
        permanent: true,
        source: "/landing-page",
      },
    ];
  },
  sassOptions: {
    functions: {
      /**
       * Copied from node_modules\@blueprintjs\core\scripts\sass-custom-functions.js
       *
       * Sass function to inline a UI icon svg and change its path color.
       *
       * Usage:
       * svg-icon("16px/icon-name.svg", (path: (fill: $color)) )
       *
       * @see https://github.com/palantir/blueprint/issues/4759#issuecomment-1112218581
       */
      "svg-icon($path, $selectors: null)": inliner(
        path.join(__dirname, "blueprint-icons/icons"), // Replaced path to icons
        {
          // minimal "uri" encoding is smaller than base64
          encodingFormat: "uri",
          // run through SVGO first
          optimize: true,
        }
      ),
    },
  },
};
