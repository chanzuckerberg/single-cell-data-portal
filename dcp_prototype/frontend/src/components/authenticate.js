import React from "react"
import {useAuth0} from "../contexts/auth0Context"
import { Box, Button, Image } from "theme-ui"

const Authenticate = () => {
  const {isAuthenticated, user, logout, loginWithRedirect} = useAuth0()

  return isAuthenticated? (
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

    </>):
    (<>
      <Button onClick={loginWithRedirect}>
        Login/Register
      </Button>
      </>)
}

export default Authenticate