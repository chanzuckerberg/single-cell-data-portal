const withImages = require("next-images");

const execSync = require("child_process").execSync;

const lastCommitCommand = "git rev-parse HEAD";

const configs = require(__dirname + "/src/configs/configs.js");
const nodeEnv = require(__dirname + "/src/common/constants/nodeEnv.js");

const { createSecureHeaders } = require("next-secure-headers");

const isProdBuild = process.env.NODE_ENV === nodeEnv.PRODUCTION;

const SCRIPT_SRC = ["'self'"];

module.exports = withImages({
  future: { webpack5: true },

  // (thuang): This is needed to ensure different Next.js servers can handle
  // the same static asset requests
  // https://github.com/vercel/next.js/issues/18389
  async generateBuildId() {
    return execSync(lastCommitCommand).toString().trim();
  },

  headers() {
    return [
      {
        headers: createSecureHeaders({
          contentSecurityPolicy: {
            directives: {
              baseUri: "'self'",
              connectSrc: [
                "'self'",
                "sentry.prod.si.czi.technology",
                configs.API_URL,
              ],
              defaultSrc: ["'self'"],
              fontSrc: ["'self'", "https://fonts.gstatic.com"],
              formAction: "'self'",
              frameAncestors: ["'none'"],
              imgSrc: ["'self'", "data:"],
              manifestSrc: ["'self'"],
              mediaSrc: ["'self'"],
              objectSrc: ["'none'"],
              reportURI:
                configs.SENTRY_DEPLOYMENT_ENVIRONMENT &&
                "https://sentry.prod.si.czi.technology/api/167/security/?sentry_key=0432f3b3ceba4bc08d28dfb61fa29707&sentry_environment=" +
                  configs.SENTRY_DEPLOYMENT_ENVIRONMENT,
              scriptSrc: isProdBuild
                ? SCRIPT_SRC
                : [...SCRIPT_SRC, "'unsafe-eval'"],
              styleSrc: [
                "'self'",
                "'unsafe-inline'",
                "https://fonts.googleapis.com",
              ],
              upgradeInsecureRequests: true,
              workerSrc: true,
            },
          },
        }),
        source: "/(.*)",
      },
    ];
  },
});
