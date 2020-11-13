// Temporarily route API calls to the backend app to http://backend
// instead of http://localhost
const configs = {
  API_URL: "http://backend:5000",
};

if (typeof module !== "undefined") module.exports = configs;
