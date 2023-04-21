export interface PayloadAction<Payload> {
  type: keyof typeof REDUCERS;
  payload: Payload;
}
export interface State {}

export const INITIAL_STATE: State = {};

export const REDUCERS = {
  dummyFunction,
};

export function reducer(state: State, action: PayloadAction<unknown>): State {
  const { type } = action;

  const handler = REDUCERS[type];

  if (!handler) {
    throw new Error(`Unknown action type: ${type}`);
  }

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  return handler(state, action as any);
}

function dummyFunction(state: State, _action: PayloadAction<null>): State {
  return state;
}
