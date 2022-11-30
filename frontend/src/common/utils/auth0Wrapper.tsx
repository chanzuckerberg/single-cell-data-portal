import React, { ReactNode } from "react";
import { useAuth0 } from "@auth0/auth0-react";

interface Props {
  children: Array<ReactNode> | ReactNode;
}

function Wrapper({ children }: Props) {
  const { isLoading, error } = useAuth0();
  if (isLoading) {
    return <div>Loading...</div>;
  }
  if (error) {
    return <div>Oops... {error.message}</div>;
  }
  return <>{children}</>;
}

export default Wrapper;
