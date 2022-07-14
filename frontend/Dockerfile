# Frontend dockerfile
FROM node:14

# install dependencies first, in a different location for easier app bind mounting for local development
# due to default /opt permissions we have to create the dir with root and change perms
RUN mkdir /opt/node_app && chown node:node /opt/node_app
WORKDIR /opt/node_app
# the official node image provides an unprivileged user as a security best practice
# but we have to manually enable it. We put it here so npm installs dependencies as the same
# user who runs the app.
# https://github.com/nodejs/docker-node/blob/master/docs/BestPractices.md#non-root-user


# -- TODO, we should try turning this back on later.
# USER node
COPY package.json package-lock.json* ./
RUN npm ci --verbose --no-optional && npm cache clean --force
ENV PATH /opt/node_app/node_modules/.bin:$PATH

RUN apt-get update && apt-get install -y \
    gconf-service gettext-base libasound2 libatk1.0-0 libatk-bridge2.0-0 libc6 libcairo2 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libgcc1 libgconf-2-4 libgdk-pixbuf2.0-0 libglib2.0-0 libgtk-3-0 libnspr4 libpango-1.0-0 libpangocairo-1.0-0 libstdc++6 libx11-6 libx11-xcb1 libxcb1 libxcomposite1 libxcursor1 libxdamage1 libxext6 libxfixes3 libxi6 libxrandr2 libxrender1 libxss1 libxtst6 ca-certificates fonts-liberation libappindicator1 libnss3 lsb-release xdg-utils wget \
    && rm -rf /var/lib/apt/lists/*

# -- TODO, we should try turning this back on later.
# ADD --chown=node . /corpora-frontend
ADD . /corpora-frontend
WORKDIR /corpora-frontend
RUN mkdir -p ./node_modules
RUN ln -sf /opt/node_app/node_modules/* /opt/node_app/node_modules/.bin ./node_modules/.
RUN if [ ! -z "$API_URL" ]; then cat src/configs/$DEPLOYMENT_STAGE.js | envsubst >src/configs/configs.js;\
    else cp src/configs/local.js src/configs/configs.js; \
    fi
RUN npm run build

ARG HAPPY_BRANCH="unknown"
ARG HAPPY_COMMIT=""
LABEL branch=${HAPPY_BRANCH}
LABEL commit=${HAPPY_COMMIT}
ENV COMMIT_SHA=${HAPPY_COMMIT}
ENV COMMIT_BRANCH=${HAPPY_BRANCH}

CMD npm run serve -- -- -p 9000
