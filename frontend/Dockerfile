# Frontend dockerfile
# Pinning to 20.10.0 to avoid GH Action lint error exit code 243
# See: https://stackoverflow.com/a/71892226/3120863
# --- Stage 1: Final Runtime Image ---
FROM node:20.10.0 AS builder

RUN npm config set registry https://registry.npmjs.org
    
# install dependencies first, in a different location for easier app bind mounting for local development
# due to default /opt permissions we have to create the dir with root and change perms
RUN mkdir /opt/node_app
WORKDIR /opt/node_app

RUN apt-get update && apt-get install -y make wget \
    gconf-service gettext-base libasound2 libatk1.0-0 libatk-bridge2.0-0 libc6 \
    libcairo2 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libgcc1 libgconf-2-4 \
    libgdk-pixbuf2.0-0 libglib2.0-0 libgtk-3-0 libnspr4 libpango-1.0-0 \
    libpangocairo-1.0-0 libstdc++6 libx11-6 libx11-xcb1 libxcb1 libxcomposite1 \
    libxcursor1 libxdamage1 libxext6 libxfixes3 libxi6 libxrandr2 libxrender1 \
    libxss1 libxtst6 ca-certificates fonts-liberation libappindicator1 libnss3 \
    lsb-release xdg-utils nano vim procps moreutils ripgrep \
    && rm -rf /var/lib/apt/lists/*

# Create the directory `make retrieve-descendants` expects before running npm ci, which runs `make retrieve-descendants`
# as a preinstall script
RUN mkdir -p src/components/common/Filter/descendant_mappings

COPY Makefile ./
COPY package*.json ./

# Run preinstall steps via `npm ci`
RUN npm ci --no-audit --prefer-offline --no-optional

# Copy app source and run build
COPY . .
ADD ./src/configs/build.js src/configs/configs.js
RUN npm run build

# --- Stage 2: Final Runtime Image ---
FROM node:20.10.0

RUN apt-get update && apt-get install -y make wget \
    gconf-service gettext-base libasound2 libatk1.0-0 libatk-bridge2.0-0 libc6 \
    libcairo2 libcups2 libdbus-1-3 libexpat1 libfontconfig1 libgcc1 libgconf-2-4 \
    libgdk-pixbuf2.0-0 libglib2.0-0 libgtk-3-0 libnspr4 libpango-1.0-0 \
    libpangocairo-1.0-0 libstdc++6 libx11-6 libx11-xcb1 libxcb1 libxcomposite1 \
    libxcursor1 libxdamage1 libxext6 libxfixes3 libxi6 libxrandr2 libxrender1 \
    libxss1 libxtst6 ca-certificates fonts-liberation libappindicator1 libnss3 \
    lsb-release xdg-utils nano vim procps moreutils ripgrep \
    && rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/node_app/node_modules/.bin:$PATH

# Copy built frontend, node_modules, and descendant mappings
COPY --from=builder /opt/node_app /opt/node_app

# Copy build output into the final working app path
WORKDIR /corpora-frontend
COPY --from=builder /opt/node_app /corpora-frontend

# If needed: config override
ADD ./src/configs/build.js src/configs/configs.js

# Link node_modules
RUN mkdir -p node_modules
RUN rm -rf node_modules && \
    cp -a /opt/node_app/node_modules ./node_modules
ARG HAPPY_BRANCH="unknown"
ARG HAPPY_COMMIT=""
LABEL branch=${HAPPY_BRANCH}
LABEL commit=${HAPPY_COMMIT}
ENV COMMIT_SHA=${HAPPY_COMMIT}
ENV COMMIT_BRANCH=${HAPPY_BRANCH}

ENTRYPOINT ["./container_init.sh"]