import React from "react"
import { handleAuthentication } from "../util/auth"

const Callback = () => {
  handleAuthentication()

  return <p>Loading...</p>
}

export default Callback
