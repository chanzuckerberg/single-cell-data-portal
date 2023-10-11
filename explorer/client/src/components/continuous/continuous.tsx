/* rc slider https://www.npmjs.com/package/rc-slider */

import React from "react";
import { connect } from "react-redux";
import HistogramBrush from "../brushableHistogram";
import Collapse from "../../util/collapse";

// @ts-expect-error ts-migrate(1238) FIXME: Unable to resolve signature of class decorator whe... Remove this comment to see the full error message
@connect((state) => ({
  // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
  schema: (state as any).annoMatrix?.schema,
}))
class Continuous extends React.PureComponent {
  // eslint-disable-next-line @typescript-eslint/explicit-module-boundary-types --- FIXME: disabled temporarily on migrate to TS.
  render() {
    /* initial value for iterator to simulate index, ranges is an object */
    // @ts-expect-error ts-migrate(2339) FIXME: Property 'schema' does not exist on type 'Readonly... Remove this comment to see the full error message
    const { schema } = this.props;
    if (!schema) return null;
    const obsIndex = schema.annotations.obs.index;
    const allContinuousNames = schema.annotations.obs.columns
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      .filter((col: any) => col.type === "int32" || col.type === "float32")
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      .filter((col: any) => col.name !== obsIndex)
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      .filter((col: any) => !col.writable) // skip user annotations - they will be treated as categorical
      // eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS.
      .map((col: any) => col.name);
    return allContinuousNames.length ? (
      <div>
        <Collapse>
          <span style={{ paddingLeft: 10 }}>Continuous</span>
          {/* eslint-disable-next-line @typescript-eslint/no-explicit-any --- FIXME: disabled temporarily on migrate to TS. */}
          {allContinuousNames.map((key: any, zebra: any) => (
            <HistogramBrush
              key={key}
              // @ts-expect-error ts-migrate(2769) FIXME: No overload matches this call.
              field={key}
              isObs
              zebra={zebra % 2 === 0}
              onGeneExpressionComplete={() => {}}
            />
          ))}
        </Collapse>
      </div>
    ) : null;
  }
}

export default Continuous;
