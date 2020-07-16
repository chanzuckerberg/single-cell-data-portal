import createAuth0Client from "@auth0/auth0-spa-js";
import Auth0Client from "@auth0/auth0-spa-js/dist/typings/Auth0Client";
import { navigate } from "gatsby";
import React, {
  createContext,
  FC,
  useContext,
  useEffect,
  useState,
} from "react";
import { User } from "../common/entities";

// (thuang): These variables are set by `env.development`, `env.production`, etc..
const config = {
  audience: process.env.BROWSER_AUDIENCE,
  // eslint-disable-next-line @typescript-eslint/camelcase
  client_id: process.env.AUTH0_CLIENT_ID as string,
  domain: process.env.AUTH0_DOMAIN as string,
  // eslint-disable-next-line @typescript-eslint/camelcase
  redirect_uri: process.env.AUTH0_CALLBACK,
};

interface Auth0Context {
  isAuthenticated: boolean;
  user: User | null;
  loading: boolean;
  handleRedirectCallback: () => void;
  loginWithRedirect: () => void;
  logout: () => void;
  getIdTokenClaims?: () => void;
  getTokenSilently?: () => Promise<string>;
}

export const Auth0Context = createContext<Auth0Context | undefined>(undefined);

export const useAuth0 = () => useContext(Auth0Context);

const Auth0Provider: FC = ({ children }) => {
  const [isAuthenticated, setIsAuthenticated] = useState(false);
  const [user, setUser] = useState<User | null>(null);
  const [auth0Client, setAuth0] = useState<Auth0Client | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    const initAuth0 = async () => {
      const auth0FromHook = await createAuth0Client(config);

      setAuth0(auth0FromHook);

      if (
        window.location.search.includes("corpora=") &&
        window.location.search.includes("state=")
      ) {
        await auth0FromHook.handleRedirectCallback();
      }

      const requestPath = window.localStorage.getItem("requestPath");
      if (requestPath) {
        window.localStorage.removeItem("requestPath");
        navigate(requestPath);
      }

      const isAuthenticated = await auth0FromHook.isAuthenticated();
      setIsAuthenticated(isAuthenticated);

      if (isAuthenticated) {
        const user = await auth0FromHook.getUser();
        console.log("user", user);
        setUser(user);
      }

      setLoading(false);
    };

    initAuth0();
  }, []);

  const loginWithRedirect = async (options?: RedirectLoginOptions) => {
    if (!auth0Client) {
      throw Error("`auth0Client` cannot be null");
    }

    if (window.localStorage.getItem("requestPath") === null) {
      window.localStorage.setItem(
        "requestPath",
        window.location.pathname + window.location.search
      );
    }

    await auth0Client?.loginWithRedirect(options);
  };

  const handleRedirectCallback = async () => {
    if (!auth0Client) {
      throw Error("`auth0Client` cannot be null");
    }

    setLoading(true);

    await auth0Client.handleRedirectCallback();

    const user = await auth0Client.getUser();

    setLoading(false);
    setIsAuthenticated(true);
    setUser(user);

    // Return to original page
    const requestPath = window.localStorage.getItem("requestPath");
    console.log("requestPath", requestPath);
    if (requestPath) {
      window.localStorage.removeItem("requestPath");
      navigate(requestPath);
    }
  };

  const logout = () => {
    if (!auth0Client) return;

    if (window.localStorage.getItem("requestPath") === null) {
      window.localStorage.setItem(
        "requestPath",
        window.location.pathname + window.location.search
      );
    }

    auth0Client.logout();
  };

  return (
    <Auth0Context.Provider
      value={{
        getIdTokenClaims: auth0Client?.getIdTokenClaims,
        getTokenSilently: auth0Client?.getTokenSilently,
        handleRedirectCallback,
        isAuthenticated,
        loading,
        loginWithRedirect,
        logout,
        user,
      }}
    >
      {children}
    </Auth0Context.Provider>
  );
};

export default Auth0Context;

export { Auth0Provider };
