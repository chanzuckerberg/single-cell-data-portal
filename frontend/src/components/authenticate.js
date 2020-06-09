import React from "react"
import { useAuth0 } from "../contexts/auth0Context"
import { Box, Flex, Image } from "theme-ui"
import StyledButton from "./styledButton"

const Authenticate = () => {
  const { isAuthenticated, loading, user, logout, loginWithRedirect } = useAuth0()

  return (
    <Box>
      {!loading && user && (
        <Flex sx={{ alignItems: "center" }}>
          <StyledButton label="Logout" onclick={logout}/>
          {!user ? (<Box>...</Box>) : (
            <Image src={user.picture}
                   width="45px"
                   height="45px"
                   padding="5px"
            />
          )}
        </Flex>
      )}
      {!loading && !isAuthenticated && (
        <>
          <StyledButton label="Login/Sign-up" onclick={loginWithRedirect}/>
        </>
      )}
    </Box>
  )

}

export default Authenticate