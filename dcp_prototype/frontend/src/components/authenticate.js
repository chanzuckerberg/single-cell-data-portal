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
      <Box>
        {!user ? (
          "..."
        ) : (
          <Image src={user.picture}
                 alt="Profile"
                 className="nav-user-profile rounded-circle"
                 width="50"  />
        )}
      </Box>
    </>):
    (<>
      <Button onClick={loginWithRedirect}>
        Login/Register
      </Button>
      <Box width="50"> </Box>
      </>)
}

export default Authenticate