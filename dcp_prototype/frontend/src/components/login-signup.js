import React from "react"
import { useAuth0 } from "../contexts/auth0Context"
import { Heading } from "theme-ui"

const LoginSignup = () => {
  const { isAuthenticated, loading, loginWithRedirect } = useAuth0()

  return (
    <>
      {!loading && !isAuthenticated && (
        <Heading as="h4" sx={{ mb: 4 }}>
          <a href="#" onClick={loginWithRedirect}>Log in</a> or <a href="#" onClick={loginWithRedirect}>sign-up</a> to
          view and download data.
        </Heading>)}
    </>
  )
}

export default LoginSignup