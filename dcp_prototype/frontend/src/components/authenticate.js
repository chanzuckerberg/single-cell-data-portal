import React from "react"
import { useAuth0 } from "../contexts/auth0Context"
import { Box, Button, Image } from "theme-ui"

const Authenticate = () => {
  const { isAuthenticated, loading, user, logout, loginWithRedirect } = useAuth0()

  return (
    <>
      {!loading && user && (
        <>
          <Button onClick={logout}>
            Logout
          </Button>

          {!user ? (<Box>...</Box>) : (
            <Image src={user.picture}
                   width="48"
                   height="48"
            />
          )}
        </>
      )}
      {!loading && !isAuthenticated && (
        <>
          <Button onClick={loginWithRedirect}>
            Login/Register
          </Button>
        </>
      )}
    </>
  )

}

export default Authenticate