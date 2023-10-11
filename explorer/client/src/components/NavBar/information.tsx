/* Core dependencies */
import React from "react";
import { connect } from "react-redux";

/* App dependencies */
import InformationMenu from "../leftSidebar/infoMenu";
import { RootState } from "../../reducers";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state: RootState) => ({
  tosURL: state.config?.parameters?.about_legal_tos,
  privacyURL: state.config?.parameters?.about_legal_privacy,
}))
class Information extends React.Component {
  handleClick = () => {
    // Dataset drawer temporarily disabled (#30).
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'dispatch' does not exist on type 'Readon... Remove this comment to see the full error message
    const { dispatch } = this.props;
    dispatch({ type: "toggle dataset drawer" });
  };

  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    const {
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'privacyURL' does not exist on type 'Read... Remove this comment to see the full error message
      privacyURL,
      // @ts-expect-error ts-migrate(2339) FIXME: Property 'tosURL' does not exist on type 'Readonly... Remove this comment to see the full error message
      tosURL,
    } = this.props;

    return (
      <div style={{ marginRight: 5, height: "100%" }}>
        <InformationMenu
          {...{
            tosURL,
            privacyURL,
          }}
        />
      </div>
    );
  }
}

export default Information;
