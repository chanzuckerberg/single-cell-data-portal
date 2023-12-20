import { isSSR } from "src/common/utils/isSSR";

export const TOO_MANY_REQUESTS_ERROR_MESSAGE_PREFIX =
  "Too many requests detected. Canceling the current and all future requests. Request count: ";

let requestCount = 0;
/**
 * (thuang): Use 5 seconds to take into account temporary spikes in requests
 */
const TIMEOUT_MS = 5 * 1000;
/**
 * (thuang): Assuming 1 request per 10ms, this will allow consecutive requests
 * for 5 seconds, which should be enough for most cases. Anything beyond that
 * is likely abnormal and needs to be investigated.
 */
const maxRequests = 500;
let timeoutExpiration = 0;
let hasReachedMaxRequests = false;

type FetchArgs = [input: RequestInfo | URL, init?: RequestInit | undefined];

export function networkGuard() {
  if (isSSR()) return;

  // Intercept network requests
  const originalFetch = window.fetch;

  window.fetch = newFetch;

  function newFetch(...args: FetchArgs) {
    requestCount++;

    if (!timeoutExpiration) {
      setTimeoutExpiration();
    }

    /**
     * (thuang): Check if the time window has expired
     * WARNING: If we've reached the max number of requests, we will no longer
     * reset the `requestCount`, and thus no new requests will be made until the
     * page is refreshed.
     */
    if (!hasReachedMaxRequests && Date.now() > timeoutExpiration) {
      resetRequestCount();
      setTimeoutExpiration();
    }

    hasReachedMaxRequests = requestCount > maxRequests;

    /**
     * (thuang): Check if the number of requests exceeds the limit
     */
    if (hasReachedMaxRequests) {
      const message = tooManyRequestsErrorMessage(requestCount);

      console.error(message);
      throw Error(message);
    }

    return originalFetch.apply(window, args);
  }
}

// Function to reset the request counter after the time window
function resetRequestCount() {
  requestCount = 0;
}

function setTimeoutExpiration() {
  timeoutExpiration = Date.now() + TIMEOUT_MS;
}

function tooManyRequestsErrorMessage(requestCount: number) {
  return TOO_MANY_REQUESTS_ERROR_MESSAGE_PREFIX + requestCount;
}
